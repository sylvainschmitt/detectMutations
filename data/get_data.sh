#!/bin/bash

from=../generateMutations/results
to=data
cp $from/base_reference/base_reference.fa $to/reference.fa
cp $from/mutated_reference/mutation.tsv $to/mutation.tsv
cp $from/base_reads1/base_reads_R2.fastq $to/base1_R2.fastq 
cp $from/base_reads1/base_reads_R1.fastq $to/base1_R1.fastq 
cp $from/base_reads2/base_reads_R1.fastq $to/base2_R1.fastq 
cp $from/base_reads2/base_reads_R2.fastq $to/base2_R2.fastq 
cp $from/mutated_reads/mixed_reads_R1.fastq $to/mutated_R1.fastq 
cp $from/mutated_reads/mixed_reads_R2.fastq $to/mutated_R2.fastq 

