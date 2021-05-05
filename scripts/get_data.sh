#!/bin/bash

cp ../generateMutations/results/base_reference/base_reference.fa results/raw_data/reference/reference.fa
cp ../generateMutations/results/mutated_reference/mutation.tsv results/raw_data/reference/mutation.tsv
cp ../generateMutations/results/base_reads1/base_reads_R2.fastq results/raw_data/reads/base1_R2.fastq 
cp ../generateMutations/results/base_reads1/base_reads_R1.fastq results/raw_data/reads/base1_R1.fastq 
cp ../generateMutations/results/base_reads2/base_reads_R1.fastq results/raw_data/reads/base2_R1.fastq 
cp ../generateMutations/results/base_reads2/base_reads_R2.fastq results/raw_data/reads/base2_R2.fastq 
cp ../generateMutations/results/mutated_reads/mixed_reads_R1.fastq results/raw_data/reads/mutated_R1.fastq 
cp ../generateMutations/results/mutated_reads/mixed_reads_R2.fastq results/raw_data/reads/mutated_R2.fastq 