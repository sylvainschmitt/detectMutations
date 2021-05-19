#!/bin/bash

from=../../generateMutations/results
mkdir reference
cp $from/reference/*.fa reference
mkdir mutations
cp $from/mutation_*/*.tsv mutations
mkdir reads
cp $from/reads_*/*.fastq reads


