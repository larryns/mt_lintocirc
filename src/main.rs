//! To align linear sequence reads to the circular mitochondrial DNA (mtDNA),
//! often what is done is to increase (double) the mitochondrial DNA reference
//! sequence. Then align the linear sequencing reads to the doubled mtDNA
//! reference. The next step--accomplished by this program--is to convert the
//! aligned reads back to an alignment of a circular mtDNA reference.

use clap::{Arg, Command};
use mt_lintocirc::convert_sam;
use noodles::sam::alignment::Record;
use noodles_util::alignment::io::reader::Builder;
use std::io;

fn main() -> io::Result<()> {
    const PROG_NAME: &str = "mt_lintocirc";
    const VERSION: &str = "0.0.1";
    const EMAIL: &str = "Larry N. Singh <larrynsingh@gmail.com>";

    let matches: clap::ArgMatches = Command::new(PROG_NAME)
            .version(VERSION)
            .author(EMAIL)
            .about("Converts BAM/SAM files mapped to doubled linear chromosome to a single linear chromosome.")
            .arg(
                Arg::new("alignmentfile")
                    .help("The input file to process")
                    .required(true) // Set to false to manually handle missing arguments
                    .index(1),
            ).get_matches();

    if let Some(filename) = matches.get_one::<String>("alignmentfile") {
        println!("Processing file: {}", filename);

        let mut reader = Builder::default().build_from_path(filename)?;

        // Process the bam file
        convert_sam::<Box<dyn Record>>(&mut reader, 16159)
    } else {
        Ok(())
    }
}
