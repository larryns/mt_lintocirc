# mt_lintocirc

## Background

The majority of software designed to align or map sequence reads to a genome were designed for the linear (nuclear) chromosomes. 
Circular chromosomes like the mitochondrial DNA (mtDNA) remain a challenge. One common approach for aligning sequence reads to the
mtDNA is to a multiple of the mtDNA reference to allow alignment of reads that map partially to the end of the reference and partially
to the beginning of the chromosome. A common approach is to double the reference genome, but I suggest using half the longest read
or insert size. 

```mt_lintocirc``` takes a BAM/SAM file aligned to an extended mtDNA reference and reverts the alignments back to the original 
linearized chromosome. In some cases, the transformation will involve splitting or breaking reads that align past the end of the linear
mtDNA. The name of the right half of the read is changed by adding the suffix *_right*, so that the read names are not duplicated. This 
program was designed with HiFi/long reads in mind, so only single end reads are handled. There is no support for paired-end reads.

## Usage

```
mt_lintocirc
    --output <output SAM file, default is to stdout>
    --alignmentfile <input alignment file>
    --ref <name of the extended mitochondrial reference; corresponds to the name in the fasta reference record>
    --reflen <length of linear reference, default is 16569>
    --targetref <output target reference name, default is chrM>
```
