//! Library (helper) modules for mt_lintocirc.

use bstr::BString;
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
        header::record::value::{map::ReferenceSequence, Map},
        io::Writer,
        Header,
    },
};
use noodles_util::alignment::io::Reader;
use std::{
    io::{self, BufRead, Write as StdWrite},
    num::NonZeroUsize,
};

// For writing we are going to use the simpler, but less efficient std::io
// writer, as opposed to the faster, but more complex tokio::async::writer.

// The return type from the convert read.
enum SplitType {
    Unchanged,           // No change to the input read
    Modified(RecordBuf), // Read modified, but not split
    Split(RecordBuf, RecordBuf),
}

/// Converts SAM records in an alignment file that were aligned to a
/// doubled circular reference genome--a reference in which the linear reference
/// genome is doubled--back to a single copy linear reference genome.
///
pub fn convert_sam<T>(
    reader: &mut Reader<Box<dyn BufRead>>,
    reflen: usize,
    bufwriter: &mut Box<dyn StdWrite>,
    refname: &BString,
    target_refname: BString,
) -> io::Result<()> {
    let mut header = reader.read_header()?;

    let reference_sequences = header.reference_sequences_mut();

    // Create an entry for chrM and add it to the end, but only if refname and target_refname are different.
    // Note that if the names don't change, we do not chagne or check the reflen.
    let mt_ref_len = NonZeroUsize::new(reflen).unwrap();
    let mt_refseq = Map::<ReferenceSequence>::new(mt_ref_len);

    // If the target_refname is the same as refname, the value is simply updated. If not,
    // a new entry is inserted and None is returned.
    if reference_sequences
        .insert(target_refname, mt_refseq)
        .is_none()
    {
        reference_sequences.swap_remove(refname);
    }

    // Open a writer to stdout. We want to lock stdout to explicitly control stdout buffering.
    let mut writer = Writer::new(bufwriter);

    // Write the header for the SAM
    writer.write_header(&header)?;

    // Loop through the SAM records.
    for result in reader.records(&header) {
        let record = result?;

        // Check if this reference is one we're interested in.

        let read_type = convert_read(&record, &header, reflen);
        match read_type {
            SplitType::Unchanged => writer.write_alignment_record(&header, &record)?,
            SplitType::Modified(read) => writer.write_alignment_record(&header, &read)?,
            SplitType::Split(left_read, right_read) => {
                // If we get a left and right read then write them separately.
                writer.write_alignment_record(&header, &left_read)?;
                writer.write_alignment_record(&header, &right_read)?;
            }
        }
    }

    // Close the writer
    writer.finish(&header)
}

fn convert_read(record: &impl Record, header: &Header, reflen: usize) -> SplitType {
    let read_name = record.name().expect("UNKNOWN read name!");

    // We are only looking for reads that are longer than `reflen`
    // First get the read alignment end.
    let Some(Ok(record_start)) = record.alignment_start() else {
        log::warn!("Record:{:?} does not have a start.", read_name);
        return SplitType::Unchanged;
    };

    // Parse the cigar vector to check if we need to split this read.
    let mut left_ref_len: usize = record_start.get();

    let name = record.name().unwrap();
    let name_str = std::str::from_utf8(name).unwrap();

    /* Often to handle a circular chromosome, the reference genome is doubled.
     * Doing so is a mistake and unnecessary because sometimes the entire
     * read will end up aligning entirely in the duplicated reference sequence.
     * To prevent this situation, the repeated part of the reference genome
     * should be no more than half of the longest expected read length. To deal
     * with this problem we subtract the reflen.
     */
    if left_ref_len >= reflen {
        log::warn!(
            "Read: {} has a start alignment: {} beyond reference.",
            name_str,
            left_ref_len
        );

        // Subtract the reflen to reset the proper alignment start.
        let mut read = RecordBuf::try_from_alignment_record(header, record).unwrap();
        let record_start = left_ref_len - reflen;
        *read.alignment_start_mut() = Position::new(record_start);

        return SplitType::Modified(read);
    }

    let mut remaining_len: usize = 0; // remaining_len is only used when there's a split
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
                    sequence_idx += curr_oper_len;
                }

                break;
            }

            // Otherwise just advance the reference counters.
            left_ref_len += oper.len();
        }

        // Move along the read?
        if oper.kind().consumes_read() {
            sequence_idx += oper.len()
        }

        // Add the operator
        left_cigar.push(oper.clone());
    }

    // Update the sequence and the quality scores for the left read.
    // Trim the sequence. In minimap2 if the alignment is a secondary alignment,
    // then there is no sequence in the secondary alignment, so nothing to split.
    let (left_sequence, right_sequence, left_quality_scores, right_quality_scores) =
        if record.sequence().len() == 0 {
            assert_eq!(
                record.quality_scores().len(),
                0,
                "Sequence for read: {} has length 0, but quality scores exist.",
                read_name
            );

            assert!(
                record.flags().unwrap().is_secondary(),
                "Sequence for read: {} has length 0, but is NOT a secondary alignment.",
                read_name
            );
            (
                Vec::<u8>::new(),
                Vec::<u8>::new(),
                Vec::<u8>::new(),
                Vec::<u8>::new(),
            )
        } else {
            // Trim the sequences and quality scores
            let mut left_seq_mut: Vec<u8> = record.sequence().iter().collect();
            let right_seq = left_seq_mut.split_off(sequence_idx);

            let mut left_qs_mut: Vec<u8> = record.quality_scores().iter().collect();
            let right_qs: Vec<u8> = left_qs_mut.split_off(sequence_idx);

            (left_seq_mut, right_seq, left_qs_mut, right_qs)
        };

    // We're done with the left cigar. Check if there was a split operator.
    // If so, add the right split part of the op to the right cigar.
    let mut right_cigar = Vec::new();

    if remaining_len > 0 {
        right_cigar.push(Op::new(curr_oper_kind, remaining_len));
    }

    // Then simply add the rest of the cigar
    while let Some(oper) = opiter.next() {
        right_cigar.push(oper.clone());
    }

    // Now create the left and right reads.
    // Create a new read alignment cloned from the record
    let mut left_read = RecordBuf::try_from_alignment_record(header, record).unwrap();

    // Update the relevant changes in the read.
    *left_read.alignment_start_mut() = Some(record_start);
    *left_read.cigar_mut() = RecordBufCigar::from(left_cigar);
    *left_read.quality_scores_mut() = RecordBufQS::from(left_quality_scores);
    *left_read.sequence_mut() = RecordBufSequence::from(left_sequence);

    // If there's nothing in the right cigar, then there's no split.
    if right_cigar.len() == 0 {
        return SplitType::Unchanged;
    }

    // Create a new read alignment cloned from the record
    let mut right_read = RecordBuf::try_from_alignment_record(header, record).unwrap();

    // Update the relevant changes in the read.
    *right_read.alignment_start_mut() = Some(Position::MIN);
    *right_read.cigar_mut() = RecordBufCigar::from(right_cigar);
    *right_read.quality_scores_mut() = RecordBufQS::from(right_quality_scores);
    *right_read.sequence_mut() = RecordBufSequence::from(right_sequence);

    // Also change the read name
    let right_name = String::from(name_str) + "_right";
    *right_read.name_mut() = Some(right_name.into());

    SplitType::Split(left_read, right_read)
}

// TESTING

#[cfg(test)]
mod tests {
    use super::*;

    use bstr::ByteSlice;
    use noodles::core::Position;
    use noodles::sam::{
        alignment::{
            record::{
                cigar::op::{Kind, Op},
                data::field::Tag,
                MappingQuality,
            },
            record_buf::{data::field::Value, RecordBuf},
        },
        header::record::value::{
            map::{Program, ReferenceSequence},
            Map,
        },
    };
    use std::{io, num::NonZeroUsize};

    const REF_LEN: usize = 1000;

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

        /* --------------------------------------------
         * | Operation | Read Index | Reference Index |
         * ------------+------------+-----------------|
         * | S20       |         20 |               0 |
         * | M30       |         50 |              30 |
         * | D5        |         50 |              35 |
         * | N5        |         50 |              40 |
         * | M90       |        140 |             130 |
         * | S10       |        150 |             130 |
         * --------------------------------------------
         *
         * Position 100 (starting at 1) is where we would like to split. At position 100
         * in the read, the reference is 10bp behind, i.e. 90bp along the reference.
         */
        let cigar: RecordBufCigar = [
            Op::new(Kind::SoftClip, 20),
            Op::new(Kind::Match, 30),
            Op::new(Kind::Deletion, 5),
            Op::new(Kind::Skip, 5),
            Op::new(Kind::Match, 90),
            Op::new(Kind::SoftClip, 10),
        ]
        .into_iter()
        .collect();

        // Generate some mapping qualities for testing
        // Convert the quality scores from string/PHRED scores to usize.
        let quality_scores = b"0123456789:;<=>?@ABCDEFGHI0123456789:;<=>?@ABCDEFGHI0123456789:;<=>?@ABCDEFGHI0123456789:;<=>?@ABCDE!FGHI0123456789:;<=>?@ABCDEFGHI0123456789:;<=>?@AB".to_vec();

        // Generate the sequence string. In this case our string is just ACGT repeated.
        // "N" marks position 100--where we are going to cut the read.
        let sequence = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA";

        // where in the reference that the read starts. See note about cigar string above.
        let record_start = REF_LEN - 90;

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
            .set_name(b"Read1".as_bstr())
            .set_template_length(100)
            .set_quality_scores(RecordBufQS::from(quality_scores))
            .set_cigar(cigar)
            .set_sequence(RecordBufSequence::from(sequence))
            .build();

        let read_type = convert_read(&sam_record, &header, REF_LEN);
        let result = match read_type {
            SplitType::Unchanged => false,
            SplitType::Modified(_) => false,
            SplitType::Split(left_read, right_read) => {
                // Check the parameters of the left and right reads.
                assert!(
                    left_read.sequence().as_ref()
                        == b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                        "left_read sequence mismatch, sequence={:?}, len={}", left_read.sequence().as_ref(), left_read.sequence().len()
                );
                assert!(
                    right_read.sequence().as_ref()
                        == b"NACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA",
                    "right_read sequence mismatch, sequence={:?}, len={}",
                    right_read.sequence().as_ref(),
                    right_read.sequence().len()
                );
                assert!(
                    left_read.quality_scores().as_ref() == b"0123456789:;<=>?@ABCDEFGHI0123456789:;<=>?@ABCDEFGHI0123456789:;<=>?@ABCDEFGHI0123456789:;<=>?@ABCDE",
                    "left_read quality mismatch, quality scores={:?}", 
                    left_read.quality_scores().as_ref()
                );
                assert!(
                    right_read.quality_scores().as_ref()
                        == b"!FGHI0123456789:;<=>?@ABCDEFGHI0123456789:;<=>?@AB",
                    "right_read quality mismatch, quality scores={:?}",
                    right_read.quality_scores().as_ref()
                );
                assert!(
                    right_read.alignment_start() == Position::new(1),
                    "right read does not start at 1"
                );
                true
            }
        };
        assert!(result, "Read not split!");

        Ok(())
    }

    #[test]
    fn test_nosplit_read() -> io::Result<()> {
        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(131072) {
            Some(n) => n,
            None => unreachable!(),
        };

        let header = noodles::sam::Header::builder()
            .set_header(Default::default())
            .add_program("nonsplit-read", Map::<Program>::default())
            .add_comment("Testing non split read function.")
            .add_reference_sequence("sq0", Map::<ReferenceSequence>::new(SQ0_LN))
            .build();

        /* --------------------------------------------
         * | Operation | Read Index | Reference Index |
         * ------------+------------+-----------------|
         * | S20       |         20 |               0 |
         * | M30       |         50 |              30 |
         * | D5        |         50 |              35 |
         * | N5        |         50 |              40 |
         * | M90       |        140 |             130 |
         * | S10       |        150 |             130 |
         * --------------------------------------------
         *
         * Position 100 (starting at 1) is where we would like to split. At position 100
         * in the read, the reference is 10bp behind, i.e. 90bp along the reference.
         */
        let cigar: RecordBufCigar = [
            Op::new(Kind::SoftClip, 20),
            Op::new(Kind::Match, 30),
            Op::new(Kind::Deletion, 5),
            Op::new(Kind::Skip, 5),
            Op::new(Kind::Match, 90),
            Op::new(Kind::SoftClip, 10),
        ]
        .into_iter()
        .collect();

        // Generate some mapping qualities for testing
        // Convert the quality scores from string/PHRED scores to usize.
        let quality_scores = b"0123456789:;<=>?@ABCDEFGHI0123456789:;<=>?@ABCDEFGHI0123456789:;<=>?@ABCDEFGHI0123456789:;<=>?@ABCDE!FGHI0123456789:;<=>?@ABCDEFGHI0123456789:;<=>?@AB".to_vec();

        // Generate the sequence string. In this case our string is just ACGT repeated.
        // "N" marks position 100--where we are going to cut the read.
        let sequence = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA";

        // where in the reference that the read starts. See note about cigar string above.
        // This starting point will make the read go over the reference by 5bp, but soft-clipping
        // doesn't consume the reference, so we won't be over the reference.
        let record_start = REF_LEN - 145;

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
            .set_name(b"Read1".as_bstr())
            .set_template_length(100)
            .set_quality_scores(RecordBufQS::from(quality_scores))
            .set_cigar(cigar)
            .set_sequence(RecordBufSequence::from(sequence))
            .build();

        let result = convert_read(&sam_record, &header, REF_LEN);
        assert!(
            matches!(result, SplitType::Unchanged),
            "Result is not none."
        );

        Ok(())
    }
}
