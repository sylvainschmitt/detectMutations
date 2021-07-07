wget https://urgi.versailles.inra.fr/download/oak/Qrob_PM1N.fa.gz
zcat Qrob_PM1N.fa.gz > Qrob_PM1N.fa
singularity pull samtools.sif oras://registry.forgemia.inra.fr/gafl/singularity/samtools/samtools:latest
./samtools.sif faidx Qrob_PM1N.fa Qrob_Chr01 -o Qrob_PM1N_Qrob_Chr01.fa
./samtools.sif faidx Qrob_PM1N.fa Qrob_Chr02 -o Qrob_PM1N_Qrob_Chr02.fa
head -n 119 Qrob_PM1N_Qrob_Chr01.fa > Qrob_PM1N_Qrob_Chr01_7k.fa
head -n 119 Qrob_PM1N_Qrob_Chr02.fa > Qrob_PM1N_Qrob_Chr02_7k.fa
cat Qrob_PM1N_Qrob_Chr01_7k.fa Qrob_PM1N_Qrob_Chr02_7k.fa > Qrob_PM1N_2Chr_7k.fa
./samtools.sif faidx Qrob_PM1N_2Chr_7k.fa
singularity pull iss.sif docker://hadrieng/insilicoseq:latest 
./iss.sif iss generate --genomes Qrob_PM1N_2Chr_7k.fa --model hiseq --n_reads 1000 --cpus 4 --o test1
./iss.sif iss generate --genomes Qrob_PM1N_2Chr_7k.fa --model hiseq --n_reads 1000 --cpus 4 --o test2
gzip *.fastq
rm Qrob_PM1N.fa Qrob_PM1N.fa.fai Qrob_PM1N_Qrob_Chr01.fa Qrob_PM1N_Qrob_Chr02.fa Qrob_PM1N_Qrob_Chr01_7k.fa Qrob_PM1N_Qrob_Chr02_7k.fa
rm samtools.sif iss.sif test1_abundance.txt test2_abundance.txt