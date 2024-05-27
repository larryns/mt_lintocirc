//! To align linear sequence reads to the circular mitochondrial DNA (mtDNA),
//! often what is done is to increase (double) the mitochondrial DNA reference
//! sequence. Then align the linear sequencing reads to the doubled mtDNA
//! reference. The next step--accomplished by this program--is to convert the
//! aligned reads back to an alignment of a circular mtDNA reference.

use bstr::BString;
use clap::{value_parser, Arg, Command};
use mt_lintocirc::convert_sam;
use noodles::sam::alignment::Record;
use noodles_util::alignment::io::reader::Builder;
use std::{
    fs::File,
    io::{self, BufWriter, Write},
};

fn main() -> io::Result<()> {
    const PROG_NAME: &str = "mt_lintocirc";
    const VERSION: &str = "0.0.1";
    const EMAIL: &str = "Larry N. Singh <larrynsingh@gmail.com>";

    let matches: clap::ArgMatches = Command::new(PROG_NAME)
            .version(VERSION)
            .author(EMAIL)
            .about("Converts BAM/SAM files mapped to doubled linear chromosome to a single linear chromosome.")
            .arg(
                Arg::new("output")
                   .short('o')
                   .long("output")
                   .required(false)
                   .help("output sam file")
            )
            .arg(
                Arg::new("alignmentfile")
                    .help("The input file to process")
                    .required(true) // Set to false to manually handle missing arguments
                    .index(1),
            ).arg(
               Arg::new("ref")
                    .short('r')
                    .long("ref")
                    .required(true)
                    .help("name of doubled mitochondrial reference") 
            ).arg(
                Arg::new("reflen")
                .short('l')
                .long("reflen")
                .required(false)
                .default_value("16569")
                .value_parser(value_parser!(u16))
                .help("length of reference, default is 16159")
            ).arg(
                Arg::new("targetref")
                .short('t')
                .long("targetref")
                .required(false)
                .default_value("chrM")
                .help("target reference sequence name")
            )
            .get_matches();

    if let Some(filename) = matches.get_one::<String>("alignmentfile") {
        log::info!("Processing file: {}", filename);

        let mut reader = Builder::default().build_from_path(filename)?;

        // Get the output file name
        let mut bufwriter: Box<dyn Write> =
            if let Some(output_filename) = matches.get_one::<String>("output") {
                let output_file = File::create_new(output_filename)?;

                Box::new(BufWriter::new(output_file))
            } else {
                Box::new(BufWriter::new(std::io::stdout().lock()))
            };

        // Get the reference name
        let refname = BString::from(matches.get_one::<String>("ref").unwrap().as_str());

        // Get the target reference name
        let target_refname =
            BString::from(matches.get_one::<String>("targetref").unwrap().as_str());

        // Get the reference name
        let target_reflen = matches.get_one::<usize>("reflen").unwrap();

        // Process the bam file
        convert_sam::<Box<dyn Record>>(
            &mut reader,
            16159,
            &mut bufwriter,
            &refname,
            target_refname,
            target_reflen,
        )
    } else {
        Ok(())
    }
}
