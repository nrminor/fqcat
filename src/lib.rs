use bio::io::fastq::{Reader, Record, Writer};
use flate2::read::GzDecoder;
use std::fs;
use std::io::{self, BufReader, BufWriter};
use std::path::PathBuf;
use zstd::stream::write::Encoder;

pub fn recode_gzip_to_zstd(input_path: &str, output_path: &str) -> io::Result<()> {
    // Create a buffer for holding 1000 records
    let mut buffer: Vec<Record> = Vec::with_capacity(1000);

    // Open the gzipped input FASTQ file
    let reader = BufReader::new(GzDecoder::new(fs::File::open(input_path)?));

    // Create FASTQ reader
    let fastq_reader = Reader::new(reader);

    // Open the Zstd compressed output FASTQ file
    let writer = BufWriter::new(Encoder::new(fs::File::create(output_path)?, 0)?);

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
                    }

                    // Clear the buffer
                    buffer.clear();
                }

                // Write any remaining records in the buffer
                for rec in &buffer {
                    fastq_writer.write_record(rec)?;
                }
            }
            Err(error) => {
                eprintln!("Encountered an error while reading: {}", error);
                continue;
            }
        }
    }

    Ok(())
}

pub struct MergePair {
    left_file: PathBuf,
    right_file: PathBuf,
}

pub struct MergeTree {
    level: usize,
    files: Vec<String>,
    merge_pairs: Option<Vec<MergePair>>,
    subtree: Option<Box<MergeTree>>,
}

pub fn build_merge_tree(
    file_list: &Vec<String>,
    level: Option<&usize>,
) -> Result<MergeTree, io::Error> {
    // handle the possibility that `level` was not specified
    let previous_level = level.unwrap_or(&0);

    // build a new tree for the provided files
    let mut tree = MergeTree {
        level: previous_level + 1,
        files: file_list.clone(),
        merge_pairs: None,
        subtree: None,
    };

    // figure out how many files will be present after merging
    let new_file_count;
    if file_list.len() % 2 == 0 {
        new_file_count = file_list.len() / 2;
    } else {
        new_file_count = (file_list.len() / 2) + 1;
    }

    // allocate a vector of that length
    let mut new_files: Vec<String> = Vec::with_capacity(new_file_count);

    // construct a list of pairs with any remainders
    let mut merge_pairs: Vec<MergePair> = Vec::new();
    let mut iter = file_list.into_iter();
    while let Some(left) = iter.next() {
        if let Some(right) = iter.next() {
            merge_pairs.push(MergePair {
                left_file: PathBuf::from(left),
                right_file: PathBuf::from(right),
            });
        } else {
            new_files.push(left.clone());
        }
    }

    // finish constructing new file list
    for (i, _) in merge_pairs.iter().enumerate() {
        let formatted_str = format!("level{}_tmp{}.fastq.zst", tree.level, i);
        new_files.push(formatted_str)
    }

    // instantiate the merge pairs and file list into the tree
    tree.merge_pairs = Some(merge_pairs);
    tree.files = new_files;

    // if there are more than two new files, construct a subtree
    if &tree.files.len() > &2 {
        tree.subtree = Some(Box::new(build_merge_tree(&tree.files, Some(&tree.level))?));
    }

    return Ok(tree);
}

fn merge_pair(pair: MergePair, output_name: String) {
    // placeholder function for the process that will handle each merge
    println!(
        "Merging {} with {} into {}",
        pair.left_file.display(),
        pair.right_file.display(),
        output_name
    );
}

pub fn process_mergepairs(pairs: Vec<MergePair>, level: usize) {
    // placeholder function for allocating a process for each merge
    println!("Merging in level {} of the merge tree", level);
    for (i, pair) in pairs.into_iter().enumerate() {
        let output_name = format!("level{}_tmp{}.fastq.zst", level, i);
        merge_pair(pair, output_name);
    }
}
