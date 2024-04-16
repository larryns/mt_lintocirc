//! To align linear sequence reads to the circular mitochondrial DNA (mtDNA),
//! often what is done is to increase (double) the mitochondrial DNA reference
//! sequence. Then align the linear sequencing reads to the doubled mtDNA
//! reference. The next step--accomplished by this program--is to convert the
//! aligned reads back to an alignment of a circular mtDNA reference.

use mt_lintocirc::convert_sam;
use noodles::sam::alignment::Record;
use noodles_util::alignment::io::reader::Builder;
use std::{env, io, path::PathBuf};

fn main() -> io::Result<()> {
    // collect turns the iterator "args" into a vector.
    let mut args = env::args();
    let filename = args
        .nth(1)
        .map(PathBuf::from)
        .expect(format!("Usage: {} <input sam/bam", args.nth(0).unwrap()).as_str());

    let mut reader = Builder::default().build_from_path(filename)?;

    // Process the bam file
    convert_sam::<Box<dyn Record>>(&mut reader, 16159)?;

    Ok(())
}
