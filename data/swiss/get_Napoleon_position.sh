# bash
cd data/swiss
zcat Napoleon_genome.fa.gz > Napoleon_genome.fa
sed -i "s/ //g" Napoleon_genome.fa
samtools faidx Napoleon_genome.fa

# R
N <- 500
mutations <- read_tsv("data/swiss/Naopleon_mutations.tsv") %>% 
   select(SNV, Contig) %>% 
  separate(Contig, c("Contig", "Pos"), convert = T) %>% 
  mutate(Start = Pos - N-1, Stop = Pos + N)
genome <- read_tsv("data/swiss/Napoleon_genome.fa.fai", 
                   col_names = c("Seq", "Size", "X1", "X2", "X3")) %>% 
  select(Seq) %>% 
  separate(Seq, paste0("X", 1:3), remove = F) %>% 
  mutate(Contig = gsub("1Quercusroburisolateech66", "", X2)) %>% 
  select(Seq, Contig)
genome %>% 
  left_join(mutations) %>% 
  filter(!is.na(SNV)) %>% 
  select(Seq, Start, Stop) %>% 
  write_tsv("data/swiss/Naopleon_mutations.bed", col_names = F)

# bash
bedtools getfasta -fi Napoleon_genome.fa -bed Naopleon_mutations.bed -fo Naopleon_mutations.fa
bwa index Qrob_PM1N.fa
makeblastdb -in Qrob_PM1N.fa -parse_seqids -blastdb_version 5 -title "Qrob_PM1N" -dbtype prot
blastn –db nt –query Naopleon_mutations –out Naopleon_mutations.out  
