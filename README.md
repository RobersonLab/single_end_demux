single_end_demux
================

Python script for demultiplexing single-end sequencing into individual FASTQ files by index.

## Use
This program expects two FASTQ formatted, phred scaled input files: a single-end Illumina sequence file and a single-end index read file. Using the index reads this program generates output FASTQ files for each index (allowing a single terminal mismatch by default). Scanning the two inputs in tandem allows the program to output the tagged reads to the correction index. Any read that does not match an index is dumped into a 'noindex' file. Any read with excessive N bases is dumped into a 'junk' file. Report of total reads, matching reads, and reads per index is printed at the end.

## Examples
    python single_end_demux seq_read.fastq index_read.fastq demux ACTCGCA,CATACCT,GCCATTG,TGACCTA
	