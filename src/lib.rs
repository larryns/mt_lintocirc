//! Library (helper) modules for mt_lintocirc.

use log::warn;
use noodles::{
    core::Position,
    sam::{
        alignment::{
            io::Write,
            record::cigar::{op::Kind, Cigar, Op},
            record_buf::{
                Cigar as RecordBufCigar, QualityScores as RecordBufQS,
                Sequence as RecordBufSequence,
            },
            Record, RecordBuf,
        },
        io::Writer,
        Header,
    },
};
use noodles_util::alignment::io::Reader;
use std::io::{self, BufRead};

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
/// convert_sam(a generic bam/sam reader, reflen);
/// ```
pub fn convert_sam<T>(reader: &mut Reader<Box<dyn BufRead>>, reflen: usize) -> io::Result<()> {
    let header = reader.read_header()?;

    // Open a writer to stdout. We want to lock stdout to explicitly control stdout buffering.
    let mut writer = Writer::new(std::io::stdout().lock());

    // Loop through the SAM records.
    for result in reader.records(&header) {
        let record = result?;

        if let Some(reads) = convert_read(&record, &header, reflen) {
            // If we get a left and right read then write them separately.
            writer.write_alignment_record(&header, &reads.0)?;
            writer.write_alignment_record(&header, &reads.1)?;
        } else {
            // Otherwise, we didn't split the read or there was a problem, in which
            // case just write the original read unchanged.
            writer.write_alignment_record(&header, &record)?;
        }
    }

    // Close the writer
    writer.finish(&header)
}

fn convert_read(
    record: &impl Record,
    header: &Header,
    reflen: usize,
) -> Option<(RecordBuf, RecordBuf)> {
    // We are only looking for reads that are longer than `reflen`
    // First get the read alignment end.
    let Some(Ok(record_start)) = record.alignment_start() else {
        let read_name = record.name().expect("UNKNOWN read name!");

        warn!("Record:{:?} does not have a start.", read_name.as_bytes());
        return None;
    };

    // Parse the cigar vector to check if we need to split this read.
    let mut left_ref_len: usize = record_start.get();
    let mut remaining_len: usize = 0;
    let mut curr_oper_kind = Kind::Skip;
    let mut left_cigar = Vec::new();
    let mut sequence_idx: usize = 0; // Used to index the sequence/qual scores

    let cigar_vec: Vec<Op> = record.cigar().iter().map(|x| x.ok().unwrap()).collect();
    let mut opiter = cigar_vec.iter();

    while let Some(oper) = opiter.next() {
        // Do we advance the reference count?
        if oper.kind().consumes_reference() {
            // If we go beyond the reference length, then we need to split this op
            if left_ref_len + oper.len() > reflen {
                // Split the oper len
                let curr_oper_len = reflen - left_ref_len;

                // Remaining len is the right size of the split from the operator.
                remaining_len = oper.len() - curr_oper_len;
                curr_oper_kind = oper.kind();

                // Add the truncated operator to the left cigar
                left_cigar.push(Op::new(oper.kind(), curr_oper_len));

                // We're done creating the left cigar. Check if we also need to advance
                // the read, then break out of the loop.
                if oper.kind().consumes_read() {
                    sequence_idx += oper.len();
                }

                break;
            }

            // Otherwise just advance the reference counters.
            left_ref_len += oper.len();

            // Add to the left cigar
            left_cigar.push(oper.clone());
        }

        // Move along the read?
        if oper.kind().consumes_read() {
            sequence_idx += oper.len()
        }
    }

    // Update the sequence and the quality scores for the left read.
    // Trim the sequence
    let mut left_sequence_mut: Vec<u8> = record.sequence().iter().collect();
    let right_sequence = left_sequence_mut.split_off(sequence_idx);

    // Trim the quality scores
    let mut left_quality_scores_mut: Vec<u8> = record.quality_scores().iter().collect();
    let right_quality_scores: Vec<u8> = left_quality_scores_mut.split_off(sequence_idx);

    // We're done with the left cigar. Check if there wa a split operator.
    // If so, add the right split part of the op to the right cigar.
    let mut right_cigar = Vec::new();

    if remaining_len > 0 {
        right_cigar.push(Op::new(curr_oper_kind, remaining_len));
    }

    // Then simply add the rest of the cigar
    while let Some(oper) = opiter.next() {
        right_cigar.push(oper.clone());
    }

    // If there's nothing in the right cigar, then there's no split.
    if right_cigar.len() == 0 {
        return None;
    }

    // Now create the left and right reads.
    // Create a new read alignment cloned from the record
    let mut left_read = RecordBuf::try_from_alignment_record(header, record).unwrap();

    // Update the relevant changes in the read.
    *left_read.alignment_start_mut() = Some(record_start);
    *left_read.cigar_mut() = RecordBufCigar::from(left_cigar);
    *left_read.quality_scores_mut() = RecordBufQS::from(left_quality_scores_mut);
    *left_read.sequence_mut() = RecordBufSequence::from(left_sequence_mut);

    // Create a new read alignment cloned from the record
    let mut right_read = RecordBuf::try_from_alignment_record(header, record).unwrap();

    // Update the relevant changes in the read.
    *right_read.alignment_start_mut() = Some(Position::MIN);
    *right_read.cigar_mut() = RecordBufCigar::from(right_cigar);
    *right_read.quality_scores_mut() = RecordBufQS::from(right_quality_scores);
    *right_read.sequence_mut() = RecordBufSequence::from(right_sequence);

    Some((left_read, right_read))
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
            .set_alignment_start(Position::new(record_start).unwrap())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::MIN)
            .set_name(Name::from(b"Read1"))
            .set_template_length(100)
            .set_quality_scores(RecordBufQS::from(quality_scores))
            .set_cigar(cigar)
            .set_sequence(RecordBufSequence::from(sequence))
            .build();

        let (left_read, right_read) = convert_read(&sam_record, &header, 50).unwrap();

        // Check the parameters of the left and right reads.

        // print results, afterwards these printlns will be changed to asserts
        println!("left = {:?}\nright = {:?}", left_read, right_read);

        assert!(false);

        Ok(())
    }
}
