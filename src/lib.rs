use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use futures::stream::FuturesUnordered;
use futures::StreamExt;
use glob::glob;
use std::fs::{remove_file, File};
use std::io::ErrorKind;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use zstd::stream::write::Encoder;
use zstd::Decoder;

pub fn find_fastqs(search_dir: &String) -> Result<Vec<String>, io::Error> {
    // Construct the search pattern
    let pattern = if search_dir.ends_with('/') {
        format!("{}*.fastq.gz", search_dir)
    } else {
        format!("{}/{}.fastq.gz", search_dir, "*")
    };

    // Initialize an empty vector to hold the paths
    let mut fastq_files: Vec<String> = Vec::new();

    // Use glob to search for files matching the pattern
    for entry in glob(&pattern)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, format!("Glob pattern error: {}", e)))?
    {
        match entry {
            Ok(path) => {
                if path
                    .file_name()
                    .and_then(|os_str| os_str.to_str())
                    .map(|file_base| !file_base.starts_with("._"))
                    .unwrap_or(false)
                {
                    fastq_files.push(path.display().to_string());
                }
            }
            Err(e) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Error recording FASTA name:\n{}", e),
                ));
            }
        }
    }

    // Check if any matching files were found
    if fastq_files.is_empty() {
        return Err(io::Error::new(ErrorKind::NotFound, "No FASTQ files found"));
    }

    Ok(fastq_files)
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
        new_file_name = if search_dir.ends_with('/') {
            format!("{}tmp{}.fastq.zst", search_dir, i)
        } else {
            format!("{}/tmp{}.fastq.zst", search_dir, i)
        };
        recode_gzip_to_zstd(&file, &new_file_name)?;
        new_files.push(new_file_name);
    }

    Ok(new_files)
}

fn recode_gzip_to_zstd(input_path: &str, output_path: &str) -> io::Result<()> {
    let file = File::open(input_path)?;

    let reader = BufReader::new(MultiGzDecoder::new(file));

    let output_file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(output_path)?;

    let mut writer = std::io::BufWriter::new(Encoder::new(output_file, 3)?.auto_finish());

    for line in reader.lines() {
        match line {
            Ok(line_content) => {
                writeln!(writer, "{}", line_content)?;
            }
            Err(e) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Error reading line:\n{}", e),
                ));
            }
        }
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
    let new_file_count: usize = if file_list.len() % 2 == 0 {
        file_list.len() / 2
    } else {
        file_list.len() / 2 + 1
    };

    // allocate a vector of that length
    let mut tmp_files: Vec<String> = Vec::with_capacity(new_file_count);

    // construct a list of pairs with any remainders
    let mut merge_pairs: Vec<MergePair> = Vec::with_capacity(new_file_count);
    let mut iter = file_list.iter();
    while let Some(left) = iter.next() {
        if let Some(right) = iter.next() {
            merge_pairs.push(MergePair {
                left_file: PathBuf::from(left),
                right_file: PathBuf::from(right),
            });
            tmp_files.push(left.to_string());
        } else {
            tmp_files.push(left.to_string());
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
    if tree.post_merge_files.len() > 1 {
        tree.subtree = Some(Box::new(build_merge_tree(
            &tree.post_merge_files,
            Some(&tree.level),
        )?));
    }

    Ok(tree)
}

async fn merge_pair(pair: MergePair) -> io::Result<()> {
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

    println!("Appending {} onto {}", file2.display(), file1.display());

    // open and buffer the file to be appended and report its extension
    let file_to_append = File::open(file2)?;
    let ext = file2.extension().expect(
        "File is incorrectly named. Be sure the file \
        contains an extension that indicates whether it is compressed.",
    );

    // decode the file to append
    if ext == "zst" {
        let decoder = zstd::Decoder::new(file_to_append)?;

        // buffer the reader, pull out the lines, and define the batch
        let read_buffer = BufReader::new(decoder);
        let mut batch: Vec<String> = Vec::with_capacity(1000);

        // Open or create the output file and create a zstd encoder for it
        let output_file = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(file1)?;
        let mut encoder = BufWriter::new(Encoder::new(output_file, 3)?.auto_finish());

        for line in read_buffer.lines() {
            let this_line = line?;
            batch.push(this_line);

            if batch.len() == 1000 {
                for batch_line in &batch {
                    writeln!(encoder, "{}", batch_line)?;
                }
                batch.clear();
            }
        }

        if !&batch.is_empty() {
            for batch_line in batch {
                writeln!(encoder, "{}", batch_line)?;
            }
        }

        println!("Merge successful; removing {:?}.", file2.display());
        if file2.display().to_string().contains("tmp") {
            remove_file(file2)?;
        }
    } else if ext == "gz" {
        let decoder = MultiGzDecoder::new(file_to_append);

        // buffer the reader, pull out the lines, and define the batch
        let read_buffer = BufReader::new(decoder);
        let mut batch: Vec<String> = Vec::with_capacity(1000);

        // Open or create the output file and create a zstd encoder for it
        let output_file = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(file1)?;
        let mut encoder = BufWriter::new(Encoder::new(output_file, 3)?.auto_finish());

        for line in read_buffer.lines() {
            let this_line = line?;
            batch.push(this_line);

            if batch.len() == 1000 {
                for batch_line in &batch {
                    writeln!(encoder, "{}", batch_line)?;
                }
                batch.clear();
            }
        }

        if !&batch.is_empty() {
            for batch_line in batch {
                writeln!(encoder, "{}", batch_line)?;
            }
        }

        println!("Merge successful; removing {:?}.", file2.display());
        if file2.display().to_string().contains("tmp") {
            remove_file(file2)?;
        }
    } else {
        // buffer the reader, pull out the lines, and define the batch
        let read_buffer = BufReader::new(file_to_append);
        let mut batch = Vec::with_capacity(1000);

        // Open or create the output file and create a zstd encoder for it
        let output_file = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(file1)?;
        let mut encoder = BufWriter::new(Encoder::new(output_file, 3)?.auto_finish());

        for line in read_buffer.lines() {
            let this_line = line?;
            batch.push(this_line);

            if batch.len() == 1000 {
                for batch_line in &batch {
                    writeln!(encoder, "{}", batch_line)?;
                }
                batch.clear();
            }
        }

        if !&batch.is_empty() {
            for batch_line in batch {
                writeln!(encoder, "{}", batch_line)?;
            }
        }

        println!("Merge successful; removing {:?}.", file2.display());
        if file2.display().to_string().contains("tmp") {
            remove_file(file2)?;
        }
    }

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
            Ok(_) => println!("Successfully appended, compressed, and cleaned {:?}", &pairs),
            Err(e) => eprintln!(
                "An error occurred when running merges in parallel:\n'{}',\nError occurred when awaiting merge completions.", e
            ),
        }
    }

    Ok(())
}

pub fn traverse_tree(tree: &MergeTree) -> io::Result<()> {
    process_mergepairs(tree.merge_pairs.clone(), tree.level)?;

    // Recur on the subtree if it exists
    if let Some(ref subtree) = tree.subtree {
        traverse_tree(subtree)?;
    }

    Ok(())
}

pub fn publish_final_fastq(readdir: &str, output_name: &str) -> io::Result<()> {
    println!("Running final conversion.");

    // handle the input
    let input_path = format!("{}/tmp0.fastq.zst", readdir);
    let open_input = File::open(&input_path)?;
    let decoder = std::io::BufReader::new(Decoder::new(open_input)?);

    // handle the output
    let open_output = File::create(output_name)?;
    let mut encoder = BufWriter::new(GzEncoder::new(open_output, Compression::default()));

    // convert to final output
    for line in decoder.lines() {
        match line {
            Ok(line_content) => {
                writeln!(encoder, "{}", line_content)?;
            }
            Err(e) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Line could not be written:\n{}", e),
                ));
            }
        }
    }

    println!(
        "Conversion complete; removing {}",
        format_args!("{}/tmp0.fastq.zst", readdir)
    );

    // remove final tmp file
    remove_file(format!("{}/tmp0.fastq.zst", readdir))?;

    Ok(())
}
