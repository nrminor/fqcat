# ReadMerger
[![Rust Build](https://github.com/nrminor/readmerger/actions/workflows/build-rust.yaml/badge.svg)](https://github.com/nrminor/readmerger/actions/workflows/build-rust.yaml) [![Open Source Files](https://github.com/nrminor/readmerger/actions/workflows/open-source-starter.yml/badge.svg)](https://github.com/nrminor/readmerger/actions/workflows/open-source-starter.yml) ![Crates.io](https://img.shields.io/crates/v/readmerger) ![Crates.io](https://img.shields.io/crates/d/readmerger)

*Never merge your reads with `find` and `cat` again!*

### Recommended Usage
Readmerger is at a very early stage of development, with many potential improvements that are detailed in [the alpha release v0.1.0](https://github.com/nrminor/readmerger/releases/tag/v0.1.0-alpha). Nevertheless, the general framework of the application is established, and it appears to be working for rapidly merging gzip- or zstd-compressed FASTQ files. We invite any interested users to use and test it with caution.

For informal use cases, we recommend users go through the following steps to install and run readmerger:

1. Make sure [the Rust toolchain](https://www.rust-lang.org/tools/install) is installed.
2. Install readmerger from [crates.io](https://crates.io/crates/readmerger) with `cargo install readmerger -v 0.1.0-alpha
`
1. Run it like so:
```
readmerger /path/to/fastq/files/ merged.fastq.gz > readmerger.log
```

This will run readmerger on your fastq files and save outputs to a log file. Those outputs include, most importantly, a pretty-printed merge tree that details the order with which files are merged.


### The Problem
Readmerger isn't finished yet, but the vision is for it to be a much faster replacement for using `cat` to merge Oxford Nanopore FASTQ files into one per barcode. Typically, Nanopore reads come out of basecalling and demultiplexing in many FASTQ files, which then must be merged. This is most often achieved with a command like this:
```
cat ~/workdir/basecall/*runid*.fastq.gz > ~/workdir/basecall/basecall.fastq.gz
```

The benefit of this method is that it uses a Unix command available on computers around the world. The trouble is, `cat` streams FASTQs into the output file one at a time, meaning that for very large sequence datasets, it can take a long time. I've had workflows where read merging is by far the slowest step. Additionally, the wildcard expansion derived from the `*` expression has an upper limit; above a certain number of files, it will error out. This leads to more verbose merging code that can look like this:
```
find . -type f -name *barcode01*.fastq.gz > fastq_list.txt
touch barcode01.fastq
for i in `cat fastq_list.txt`;
do
    zcat $i >> barcode01.fastq
done
gzip barcode01.fastq
```

While this code solves the wildcard expansion problem, it still merges one FASTQ at a time, and to avoid compression corruption, it has to decompress and recompress as well. To solve the wildcard expansion issue, we've ended up with code that's even slower!

### Readmerger's approach
The vision for `readmerger` is that it will solve these problems and be as easy to run as `cat` (while also being easy to install). The tool will use what I call a "merge tree", or "hierarchical merging." As an example, say there are eight FASTQs in a barcode directory. Readmerger will take each of the four pairs of FASTQs and merge each pair in parallel, resulting in 4 merged FASTQs. Then, it will do the same thing, merging each pair in parallel and outputting two FASTQs, at which point it will merge the last two and give you your final FASTQ. Like a March Madness Bracket, readmerger will start with the many tips of the tree of FASTQs and end up with one winner.

Importantly, the tool will avoid redundant file reading and writing by using an append model: in the first "round" of the merge tree, it will transfer one of the FASTQs reads into a new FASTQ and then append reads from other FASTQs to that file for the rest of the merging process.
