//! Library (helper) modules for mt_lintocirc.

use log::warn;
use noodles::bam;
use noodles::core::Position;
use noodles::sam;
use noodles::sam::alignment::record::{cigar::Op, Data};
use noodles::sam::alignment::record::{Cigar, Name, Record};
use noodles::sam::alignment::record_buf::{Cigar as CigarBuf, RecordBuf};
use noodles::sam::Header;
use noodles::util::alignment;
use std::{i32, io};

// TODO: Look up alignment_span in noodles for how to parse the cigar string.

// For writing we are going to use the simpler, but less efficient std::io
// writer, as opposed to the faster, but more complex tokio::async::writer.

/// Converts SAM records in an alignment file that were aligned to a
/// doubled circular reference genome--a reference in which the linear reference
/// genome is doubled--back to a single copy linear reference genome.
///
/// # Example
///
/// ```
/// use mt_lintocirc::convert_bam;
/// convert_bam(filename, reflen);
/// ```
pub fn convert_bam(filename: &str, reflen: usize) -> io::Result<()> {
    const MAX_CIGAR_LEN: usize = 1 << 16;

    // Open a reader to the file given.
    let mut reader = bam::io::reader::Builder.build_from_path(filename)?;
    let header = reader.read_header()?;

    // Open a writer to stdout. We want to lock stdout to explicitly control stdout buffering.
    let mut writer = alignment::io::writer::builder::Builder::default()
        .set_format(alignment::io::Format::Sam)
        .build_from_writer(std::io::stdout().lock)?;

    // Loop through the SAM records.
    for result in reader.records() {
        let record = result?;

        // We are only looking for reads that are longer than `reflen`
        // First get the read alignment end.
        let Some(Ok(record_start)) = record.alignment_start() else {
            let read_name = record.name().expect("UNKNOWN read name!");

            warn!("Record:{:?} does not have a start.", read_name.as_bytes());
            continue;
        };

        // We need to build two new reads. Based on how this library is designed,
        // we can't just modify one read. First, we build the left read. We will
        // need a new (shorter) sequence, qualities, cigar string, template length,
        let mut left_ref_len: usize = record_start.get();
        let mut left_read_len: usize = 0;

        // Traverse the read to check if we go past the end of the reference genome.
        for opiter in record.cigar().iter() {
            let oper = opiter?;

            if oper.kind().consumes_reference() {
                left_ref_len += 1;

                // check if we are at the boundary of the reference
                if left_ref_len > reflen {
                    // We need to split the read.
                    break;
                }
            }

            if oper.kind().consumes_read() {
                left_read_len += 1;
            }
        }

        // At this point, we either need to split the read or just write it out
        if left_ref_len > reflen {
            let (left_read, right_read) = split_read(&record, &header, left_read_len, record_start);
            writer.write_record(&header, &left_read)?;
            writer.write_record(&header, &right_read)?;
        } else {
            writer.write_record(&header, &record)?;
        }
    }

    // Close the writer
    writer.finish(&header)?;

    Ok(())
}

// Split a read into two parts, and return both parts.
fn split_read(
    record: &noodles::bam::Record,
    header: &Header,
    at: usize,
    record_start: Position,
) -> (RecordBuf, RecordBuf) {
    // Split the cigar string for the two reads
    let left_cigar: CigarBuf = record
        .cigar()
        .iter()
        .take(at)
        .map(|x| x.ok().unwrap())
        .collect();
    let right_cigar = record.cigar().iter().skip(at).collect();

    // Trim the quality scores
    let left_quality_scores = record.quality_scores().into().iter().take(at);
    let left_sequence = record.sequence().into().iter().take(at);

    // Create a new read alignment cloned from the record
    let mut left_read = RecordBuf::try_from_alignment_record(header, record).unwrap();

    // Update the relevant changes in the read.
    *left_read.alignment_start_mut() = Some(record_start);
    *left_read.cigar_mut() = left_cigar;
    *left_read.quality_scores_mut() = left_quality_scores;
    *left_read.sequence_mut() = left_sequence;
    *left_read.template_length_mut() = i32::try_from(at).unwrap();

    // Now the right read
    let right_cigar = record
        .cigar()
        .iter()
        .skip(at)
        .map(|x| x.ok().unwrap())
        .collect();

    // Trim the quality scores
    let left_quality_scores = record.quality_scores().into().iter().take(at);
    let left_sequence = record.sequence().into().iter().take(at);

    // Create a new read alignment cloned from the record
    let mut right_read = RecordBuf::try_from_alignment_record(header, record).unwrap();

    let right_quality_scores = record.quality_scores().into().iter().skip(at);
    let right_sequence = record.sequence().into().iter().skip(at);

    // Update the relevant changes in the read.
    *right_read.alignment_start_mut() = Some(Position::MIN);
    *right_read.cigar_mut() = right_cigar;
    *right_read.quality_scores_mut() = right_quality_scores;
    *right_read.sequence_mut() = right_sequence;
    *right_read.template_length_mut() = right_read.template_length() - i32::try_from(at).unwrap();

    (left_read, right_read)
}

// TESTING

#[cfg(test)]
mod tests {
    use super::*;

    const MTDNA_REF_LEN: usize = 16159;

    #[test]
    fn test_split_cigar() {
        let cigar = "20S30M5D5N40M10S"
            .parse::<cigar::Cigar>()
            .expect("FATAL: test_split_failed due to Cigar parse failure.");

        let (left, right) = split_cigar(&cigar, MTDNA_REF_LEN - 10, MTDNA_REF_LEN);

        println!("left: {:?}", left);
        println!("right: {:?}", right);

        let left_cigar: cigar::Cigar = left.into_iter().collect::<cigar::Cigar>();
        let right_cigar: cigar::Cigar = right.into_iter().collect::<cigar::Cigar>();

        assert_eq!("20S10M".parse::<cigar::Cigar>(), Ok(left_cigar));
        assert_eq!("20M5D5N40M10S".parse::<cigar::Cigar>(), Ok(right_cigar));
    }
}
