use clap::Parser;
use glob::glob;
use readmerger::process_fastq_gzip_to_zstd;
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
    let output_path = args.output_name;

    // Construct the search pattern
    let pattern = format!("{}/*.fastq.gz", input_dir);

    // Initialize an empty vector to hold the paths
    let mut fastq_files: Vec<String> = Vec::new();

    // Use glob to search for files matching the pattern
    for entry in
        glob(&pattern).expect("Failed to find any FASTQ files in the provided directory pattern")
    {
        match entry {
            Ok(path) => {
                let path_str = path.display().to_string();
                fastq_files.push(path_str);
            }
            Err(e) => println!("{:?}", e),
        }
    }

    // loop through each FASTQ and convert them to ZSTD
    for fastq in fastq_files {
        process_fastq_gzip_to_zstd(&fastq, &output_path)?;
    }

    println!("Processing completed successfully.");
    Ok(())
}
