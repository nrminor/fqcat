use bio::io::fastq::{Reader, Record, Writer};
use flate2::read::GzDecoder;
use std::fs;
use std::io::{self, BufReader, BufWriter};
use zstd::stream::write::Encoder;

pub fn process_fastq_gzip_to_zstd(input_path: &str, output_path: &str) -> io::Result<()> {
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
