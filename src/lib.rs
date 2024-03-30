//! Library (helper) modules for mt_lintocirc.

use log::warn;
use noodles::bam;
use noodles::core::Position;
use noodles::sam;
use noodles::sam::alignment::record::cigar::Op;
use noodles::sam::alignment::record::{Cigar, Name, Record};
use noodles::sam::alignment::record_buf::RecordBuf;

use std::io;

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
    let stdout = io::stdout().lock();
    let mut writer = sam::io::Writer::new(stdout);

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
            let (left_read, right_read) = split_read(left_read_len);
            writer.write_alignment_record(&header, &left_read)?;
            writer.write_alignment_record(&header, &right_read)?;
        } else {
            writer.write_alignment_record(&header, &record)?;
        }

        if alignment_span >= reflen {
            // Traverse the cigar operations and split the cigar
            // string into two.

            let (left, right) = split_cigar(&record.cigar(), record_start, reflen);

            // Now update the record with cigar1, and then clone
            // the record and update the alignment length and cigar with
            // cigar2. The simplest way as far as I can tell to change the
            // record cigar, unfortunately, is to create a completely new record.
            let qs_left = &record.quality_scores().as_ref()[0..left.len()];
            let seq_left = &record.sequence().as_ref()[0..left.len()];
            let cigar_left = Cigar::from(left);

            let left_rec = RecordBuf::builder()
                .set_alignment_start(record_start)
                .set_cigar(left.into_iter().collect())
                .set_data(record.data().clone())
                .set_flags(record.flags())
                .set_mapping_quality(record.mapping_quality().unwrap())
                .set_quality_scores(qs_left)
                .set_read_name(*record.read_name().unwrap())
                .set_reference_sequence_id(record.reference_sequence_id().unwrap())
                .set_sequence(seq_left)
                .set_template_length(left.len() as i32)
                .build();

            // Repeat for the right side
            let readname_str: &str = record.read_name().unwrap().as_ref();
            let readname_right = format!("{readname_str}_right");
            let new_readname = Name::from_str(readname_right.as_str()).unwrap();

            /*
            let right_rec = sam::alignment::Record::builder()
                .set_alignment_start(Position::new(1).unwrap())
                .set_cigar(right.into_iter().collect())
                .set_data(record.data().clone())
                .set_flags(record.flags())
                .set_mapping_quality(record.mapping_quality().unwrap())
                .set_quality_scores(sam::record::QualityScores::from(new_qs))
                .set_read_name(new_readname)
                .set_reference_sequence_id(record.reference_sequence_id().unwrap())
                .set_sequence(noodles_sam::record::Sequence::from(*new_sequence))
                .set_template_length(right.len() as i32)
                .build();
            */
        } else {
            // Simply output the unchanged record.
            writer.write_alignment_record(&header, &record)?;
        }
    }

    // Close the writer
    writer.finish(&header)?;

    Ok(())
}

/// Splits a SAM CIGAR record aligned across the boundary of a doubled reference
/// genome back into a single copy linear genome. This process involves breaking
/// the CIGAR string right at the boundary and adjusting the CIGAR Op at the
/// breakpoint.
fn split_cigar(cigar: &dyn Cigar, align_start: Position, reflen: usize) -> (Vec<Op>, Vec<Op>) {
    // Loop through the CIGAR operations until we get to the operation
    // that crosses the boundary of the reference. As we traverse `cigar`
    // build a copy for cigar1.
    let mut currlen = align_start; // current length of reference
    let mut left_op_vec: Vec<Op> = vec![];

    // First process the left side. I couldn't find a way to add an
    // op to an existing vec, so we will have to destruct and rebuild
    // the cigar strings. This approach seems to be very inefficient,
    // but I don't see another way.
    let mut op_iter = cigar.iter();
    let mut target_op: Option<&Op> = None;

    while let Some(op) = op_iter.next() {
        let oplen = op
            .kind()
            .consumes_reference()
            .then_some(op.len())
            .unwrap_or(0);

        // Are we at the Op that crosses the boundary?
        if currlen + oplen >= reflen {
            target_op = Some(op);
            break;
        }

        // update the running reference length
        currlen = currlen + oplen;

        // Add the operations to the left side.
        left_op_vec.push(*op);
    }

    // target_op better not be None or we have a logic error.
    // split_op is the op we need to split.
    let split_op = target_op.expect("Logic error!");

    // Compute the new left "op"
    let left_op: cigar::Op = cigar::Op::new(split_op.kind(), reflen - currlen);
    left_op_vec.push(left_op);

    // Now the right op
    let right_op: cigar::Op = cigar::Op::new(split_op.kind(), split_op.len() - left_op.len());
    let mut right_op_vec: Vec<cigar::Op> = vec![right_op];

    // Add the rest of the ops.
    while let Some(op) = op_iter.next() {
        right_op_vec.push(*op);
    }

    // Finally convert the vector of Ops to a cigar string
    (left_op_vec, right_op_vec)
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
