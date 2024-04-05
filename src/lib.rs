//! Library (helper) modules for mt_lintocirc.

use log::warn;
use noodles::bam;
use noodles::core::Position;
use noodles::sam::alignment::io::Write;
use noodles::sam::alignment::record_buf::{
    Cigar as RecordBufCigar, QualityScores as RecordBufQS, Sequence as RecordBufSequence,
};
use noodles::sam::alignment::RecordBuf;
use noodles::sam::io::Writer;
use noodles::sam::Header;
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
/// use mt_lintocirc::convert_sam;
/// convert_sam(filename, reflen);
/// ```
pub fn convert_sam(filename: &str, reflen: usize) -> io::Result<()> {
    // Open a reader to the file given.
    let mut reader = bam::io::reader::Builder.build_from_path(filename)?;
    let header = reader.read_header()?;

    // Open a writer to stdout. We want to lock stdout to explicitly control stdout buffering.
    let mut writer = Writer::new(std::io::stdout().lock());

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
            writer.write_alignment_record(&header, &left_read)?;
            writer.write_alignment_record(&header, &right_read)?;
        } else {
            writer.write_alignment_record(&header, &record)?;
        }
    }

    // Close the writer
    writer.finish(&header)?;

    Ok(())
}

// Split a read into two parts, and return both parts.
fn split_read(
    record: &impl noodles::sam::alignment::Record,
    header: &Header,
    at: usize,
    record_start: Position,
) -> (RecordBuf, RecordBuf) {
    // Split the cigar string for the two reads
    let left_cigar: RecordBufCigar = record
        .cigar()
        .iter()
        .take(at)
        .map(|x| x.ok().unwrap())
        .collect();
    let right_cigar: RecordBufCigar = record.cigar().iter().skip(at).map(|x| x.unwrap()).collect();

    // Trim the quality scores
    let mut left_quality_scores_mut: Vec<u8> = record.quality_scores().iter().collect();
    let right_quality_scores: Vec<u8> = left_quality_scores_mut.split_off(at);

    // Trim the sequence
    let mut left_sequence_mut: Vec<u8> = record.sequence().iter().collect();
    let right_sequence = left_sequence_mut.split_off(at);

    // Create a new read alignment cloned from the record
    let mut left_read = RecordBuf::try_from_alignment_record(header, record).unwrap();

    // Update the relevant changes in the read.
    *left_read.alignment_start_mut() = Some(record_start);
    *left_read.cigar_mut() = left_cigar;
    *left_read.quality_scores_mut() = RecordBufQS::from(left_quality_scores_mut);
    *left_read.sequence_mut() = RecordBufSequence::from(left_sequence_mut);
    *left_read.template_length_mut() = i32::try_from(at).unwrap();

    // Create a new read alignment cloned from the record
    let mut right_read = RecordBuf::try_from_alignment_record(header, record).unwrap();

    // Update the relevant changes in the read.
    *right_read.alignment_start_mut() = Some(Position::MIN);
    *right_read.cigar_mut() = right_cigar;
    *right_read.quality_scores_mut() = RecordBufQS::from(right_quality_scores);
    *right_read.sequence_mut() = RecordBufSequence::from(right_sequence);
    *right_read.template_length_mut() = right_read.template_length() - i32::try_from(at).unwrap();

    (left_read, right_read)
}

// TESTING

#[cfg(test)]
mod tests {
    use super::*;

    use noodles::core::Position;
    use noodles::sam::{
        alignment::{
            record::{
                cigar::op::{Kind, Op},
                data::field::Tag,
                MappingQuality,
            },
            record_buf::{data::field::Value, Name, RecordBuf},
        },
        header::record::value::{
            map::{Program, ReferenceSequence},
            Map,
        },
    };
    use std::{io, num::NonZeroUsize};

    const REF_LEN: usize = 1000;
    const READ_LEN: usize = 100;

    #[test]
    fn test_split_read() -> io::Result<()> {
        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(131072) {
            Some(n) => n,
            None => unreachable!(),
        };

        let header = noodles::sam::Header::builder()
            .set_header(Default::default())
            .add_program("split-read", Map::<Program>::default())
            .add_comment("Testing split_read function.")
            .add_reference_sequence("sq0", Map::<ReferenceSequence>::new(SQ0_LN))
            .build();

        let cigar: RecordBufCigar = [
            Op::new(Kind::SoftClip, 20),
            Op::new(Kind::Match, 30),
            Op::new(Kind::Deletion, 5),
            Op::new(Kind::Skip, 5),
            Op::new(Kind::Match, 40),
            Op::new(Kind::SoftClip, 10),
        ]
        .into_iter()
        .collect();

        // Generate some mapping qualities for testing
        let quality_scores: Vec<_> = (0..READ_LEN)
            .map(|i| (((i as f32) / 10.0 * 3.0 + 10.0) as u8))
            .collect();

        // Generate the sequence string
        let sequence = b"ACGT".repeat(READ_LEN / 4);

        // where in the reference, that the read starts
        let record_start = REF_LEN - 30;

        // Create a SAM record to split
        let sam_record = RecordBuf::builder()
            .set_data(
                [
                    (Tag::READ_GROUP, Value::from("rg0")),
                    (Tag::ALIGNMENT_HIT_COUNT, Value::UInt8(1)),
                ]
                .into_iter()
                .collect(),
            )
            .set_alignment_start(Position::new(REF_LEN - 30).unwrap())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::MIN)
            .set_name(Name::from(b"Read1"))
            .set_template_length(100)
            .set_quality_scores(RecordBufQS::from(quality_scores))
            .set_cigar(cigar)
            .set_sequence(RecordBufSequence::from(sequence))
            .build();

        let (left_read, right_read) = split_read(
            &sam_record,
            &header,
            50,
            Position::new(record_start).unwrap(),
        );

        // print results, afterwards these printlns will be changed to asserts
        println!("left = {:?}\nright = {:?}", left_read, right_read);

        Ok(())
    }
}
