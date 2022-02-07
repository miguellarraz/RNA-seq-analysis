## RNA-seq-analysis

This code is an example of a pipeline that can be used to calculate fold-change data from a RNA-seq experiment, using exclusively bash.

Initially, an index file is produced for the reference genome using bowtie, and fastqc is run to produce a report showing the quality of each read.
Each read is aligned, flagging poorly aligned reads. Bam files are produced, which are sorted and indexed, and the coverage is calculated with bedtools (this is done to find the counts in each sample to produce fold-change data).
Finally, a file is produce representing the fold change between all samples.
