//! To align linear sequence reads to the circular mitochondrial DNA (mtDNA),
//! often what is done is to increase (double) the mitochondrial DNA reference
//! sequence. Then align the linear sequencing reads to the doubled mtDNA
//! reference. The next step--accomplished by this program--is to convert the
//! aligned reads back to an alignment of a circular mtDNA reference.

use mt_lintocirc::convert_sam;
use std::{env, io};

fn main() -> io::Result<()> {
    // collect turns the iterator "args" into a vector.
    let args: Vec<String> = env::args().collect();

    if args.len() != 2 {
        eprintln!("Usage: {} <input bam>.", args[0]);
        return Ok(());
    }

    // Process the bam file
    convert_sam(&args[1], 16159)?;

    Ok(())
}
