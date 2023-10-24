use clap::Parser;
use readmerger::{build_merge_tree, find_fastqs, prepare_for_merges, traverse_tree};
use std::io;

// Search for a pattern in a file and display the lines that contain it.
#[derive(Parser)]
struct Cli {
    #[clap(value_enum)]
    readdir: String,
    output_name: String,
}

fn main() -> io::Result<()> {
    let args = Cli::parse();
    let input_dir = args.readdir;
    let _output_path = args.output_name;

    // find input FASTQ files in the provided directory
    let fastq_files = find_fastqs(&input_dir)?;

    // determine which files will be appended to while re-compressing with Zstandard
    let prepped_files = prepare_for_merges(fastq_files)?;

    // construct merge tree
    let merge_tree = build_merge_tree(&prepped_files, None)?;

    // traverse the tree and merge file pairs until none remain
    _ = traverse_tree(&merge_tree)?;

    Ok(())
}
