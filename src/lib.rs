use bio::io::fastq::{Reader, Record, Writer};
use flate2::read::GzDecoder;
use futures::stream::FuturesUnordered;
use futures::StreamExt;
use glob::glob;
use std::fs;
use std::io::ErrorKind;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::PathBuf;
use zstd::stream::write::Encoder;

/*
NOTES TO SELF:

Will need to amend the functions below to retain the tmp fastq name from function
prepare_for_merges() throughout the merge tree and always append to it. Will also need
to code in behavior for if only two fastqs remain; in that case, merge them, convert
them to gzip, and apply the currently unused output filename.

*/

pub fn find_fastqs(search_dir: &String) -> Result<Vec<String>, io::Error> {
    // Construct the search pattern
    let pattern = format!("{}/*.fastq.gz", search_dir);

    // Initialize an empty vector to hold the paths
    let mut fastq_files: Vec<String> = Vec::new();

    // Use glob to search for files matching the pattern
    for entry in
        glob(&pattern).expect("Failed to find any FASTQ files in the provided directory pattern")
    {
        match entry {
            Ok(path) => {
                let path_str = path.display().to_string();
                let file_base = path.file_name().unwrap().to_str().unwrap();
                if !file_base.starts_with("._") {
                    fastq_files.push(path_str);
                }
            }
            Err(e) => println!("{:?}", e),
        }
    }

    // Check if any matching files were found
    if fastq_files.is_empty() {
        return Err(io::Error::new(ErrorKind::NotFound, "No FASTQ files found"));
    }

    return Ok(fastq_files);
}

pub fn prepare_for_merges(
    file_list: Vec<String>,
    search_dir: &String,
) -> Result<Vec<String>, io::Error> {
    let mut new_files = Vec::with_capacity(file_list.len());
    let mut new_file_name: String;

    for (i, file) in file_list.into_iter().enumerate() {
        if i % 2 != 0 {
            new_files.push(file);
            continue;
        }
        new_file_name = format!("{}/tmp{}.fastq.zst", search_dir, i);
        recode_gzip_to_zstd(&file, &new_file_name)?;
        new_files.push(new_file_name);
    }

    return Ok(new_files);
}

fn recode_gzip_to_zstd(input_path: &str, output_path: &str) -> io::Result<()> {
    // Create a buffer for holding 1000 records
    let mut buffer: Vec<Record> = Vec::with_capacity(1000);

    // Open the gzipped input FASTQ file
    let reader = BufReader::new(GzDecoder::new(fs::File::open(input_path)?));

    // Create FASTQ reader
    let fastq_reader = Reader::new(reader);

    // Open the Zstd compressed output FASTQ file
    let writer = BufWriter::new(Encoder::new(fs::File::create(output_path)?, 3)?);

    // Create FASTQ writer
    let mut fastq_writer = Writer::new(writer);

    // Loop over records
    for result in fastq_reader.records() {
        match result {
            Ok(record) => {
                buffer.push(record);

                if buffer.len() == 1000 {
                    // Write the buffered records to the Zstd compressed output file
                    for rec in &buffer {
                        fastq_writer.write_record(rec)?;
                        fastq_writer.flush()?;
                    }

                    // Clear the buffer
                    buffer.clear();
                }
            }
            Err(error) => {
                eprintln!("Encountered an error while reading: {}", error);
                continue;
            }
        }
    }

    // Write any remaining records in the buffer
    for rec in &buffer {
        fastq_writer.write_record(rec)?;
        fastq_writer.flush()?;
    }

    Ok(())
}

#[derive(Clone, Debug)]
pub struct MergePair {
    left_file: PathBuf,
    right_file: PathBuf,
}

#[derive(Debug)]
pub struct MergeTree {
    level: usize,
    post_merge_files: Vec<String>,
    merge_pairs: Vec<MergePair>,
    subtree: Option<Box<MergeTree>>,
}

pub fn build_merge_tree(
    file_list: &Vec<String>,
    level: Option<&usize>,
) -> Result<MergeTree, io::Error> {
    // handle the possibility that `level` was not specified
    let previous_level = level.unwrap_or(&0);

    // figure out how many files will be present after merging
    let new_file_count: usize;
    if file_list.len() % 2 == 0 {
        new_file_count = file_list.len() / 2;
    } else {
        new_file_count = (file_list.len() / 2) + 1;
    }

    // allocate a vector of that length
    let mut tmp_files: Vec<String> = Vec::with_capacity(new_file_count);

    // construct a list of pairs with any remainders
    let mut merge_pairs: Vec<MergePair> = Vec::with_capacity(new_file_count);
    let mut iter = file_list.into_iter();
    while let Some(left) = iter.next() {
        if let Some(right) = iter.next() {
            merge_pairs.push(MergePair {
                left_file: PathBuf::from(left),
                right_file: PathBuf::from(right),
            });
            tmp_files.push(left.clone());
        } else {
            tmp_files.push(left.clone());
        }
    }

    // build a new tree for the provided files
    let mut tree = MergeTree {
        level: previous_level + 1,
        post_merge_files: tmp_files,
        merge_pairs,
        subtree: None,
    };

    // if there are more than two new files, construct a subtree
    if &tree.post_merge_files.len() > &1 {
        tree.subtree = Some(Box::new(build_merge_tree(
            &tree.post_merge_files,
            Some(&tree.level),
        )?));
    }

    return Ok(tree);
}

async fn merge_pair(pair: MergePair) -> io::Result<()> {
    // placeholder function for the process that will handle each merge

    // convert the left file in the pair to a string
    let left_file = pair
        .left_file
        .to_str()
        .ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, "Left file was None"))?;

    // based on whether the left file is an intermediate file that has
    // resulted from a previous merge, bind file paths so that the correct
    // file is appended to
    let (file1, file2) = if left_file.contains("tmp") {
        (&pair.left_file, &pair.right_file)
    } else {
        (&pair.right_file, &pair.left_file)
    };

    // open and buffer the file to be appended
    let file_to_append = fs::File::open(&file2)?;

    // decode the file to append based on whether or how it is compressed
    let reader: Box<dyn Read> = if file2.ends_with(".zst") {
        Box::new(zstd::Decoder::new(file_to_append)?)
    } else if file2.ends_with(".gz") {
        Box::new(GzDecoder::new(file_to_append))
    } else {
        Box::new(file_to_append)
    };

    // buffer the reader and define the batch
    let read_buffer = BufReader::new(reader);
    let mut batch: Vec<String> = Vec::new();

    // Open or create the output file and create a zstd encoder for it
    let output_file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(file2)?;
    let mut encoder = Encoder::new(output_file, 0)?;

    for line in read_buffer.lines() {
        let line = line?;
        batch.push(line);

        if batch.len() >= 1000 {
            for line in &batch {
                writeln!(encoder, "{}", line)?;
            }
            batch.clear();
        }
    }

    if !batch.is_empty() {
        for line in batch {
            writeln!(encoder, "{}", line)?;
        }
    }

    encoder.finish()?;

    println!("Merging {} with {}", file1.display(), file2.display());

    Ok(())
}

#[tokio::main]
async fn process_mergepairs(pairs: Vec<MergePair>, level: usize) -> io::Result<()> {
    println!("Processing file pairs at level {} of the merge tree", level);

    let mut futures = FuturesUnordered::new();

    for pair in &pairs {
        futures.push(merge_pair(pair.clone()));
    }

    while let Some(result) = futures.next().await {
        match result {
            Ok(_) => println!("Successfully appended and compressed."),
            Err(e) => eprintln!("An error occurred: {}", e),
        }
    }

    Ok(())
}

pub fn traverse_tree(tree: &MergeTree, output_name: &str) -> io::Result<()> {
    process_mergepairs(tree.merge_pairs.clone(), tree.level.clone())?;

    // Recur on the subtree if it exists
    if let Some(ref subtree) = tree.subtree {
        traverse_tree(subtree, &output_name)?;
    }

    Ok(())
}
