library(tidyverse)
library(vroom)
library(Biostrings)
library(genomeIntervals)

mutations <- vroom(snakemake@input[[1]])
te_gff <- readGff3(snakemake@input[[2]], quiet = T)
te_tab <- vroom(snakemake@input[[2]], skip = 1,
                   col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))
ref <- snakemake@config[["reference"]]
out <-  snakemake@output[[1]]

# mutations <- vroom("results/mutations/mutations_filtered.tsv")
# te_gff <- readGff3("data/genomes/FsylCur4_refTEs_wclassif_wreliable.gff", quiet = T)
# te_tab <- vroom("data/genomes/FsylCur4_refTEs_wclassif_wreliable.gff", skip = 1,
#                    col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))
# ref <- "Fagus_sylvatica_v3"

mutations <- mutations %>% 
  filter(reference == ref) %>% 
  mutate(SNV = paste0(CHROM, "#", POS))

mutationsI <- new(
  "Genome_intervals_stranded",
  as.matrix(mutations[c("POS", "POS")]),
  closed = T,
  annotation = data.frame(
    seq_name = mutations$CHROM,
    inter_base = F,
    strand = factor("+", levels=c("+", "-") )
  )
)

overlaps <- interval_overlap(mutationsI, te_gff)

te <- lapply(1:length(overlaps), function(i){
  if(length(overlaps[[i]]) > 0)
    bind_cols(
      mutations[i,],
      te_tab[overlaps[[i]],] %>% 
        dplyr::select(-seqid, -source) %>% 
        summarise(type = paste(type, collapse = ";"), start = min(start), end = max(end), score = dplyr::first(score),
                  strand = dplyr::first(strand), phase = dplyr::first(phase), attributes = paste(attributes, collapse = ";")))
  
}) %>% bind_rows() %>% 
  dplyr::select(SNV, type, start, end, score, strand, phase, attributes) %>% 
  dplyr::rename(te_type = type, te_start = start, te_end = end, te_score = score, te_strand = strand,
                te_phase = phase, te_attributes = attributes)

write_tsv(te, out)

