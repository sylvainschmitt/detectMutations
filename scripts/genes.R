library(tidyverse)
library(vroom)
library(Biostrings)
library(genomeIntervals)

mutations <- vroom(snakemake@input[[1]])
genes_gff <- readGff3(snakemake@input[[2]], quiet = T)
genes_tab <- vroom(snakemake@input[[2]], skip = 1,
                   col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))
annotation <- vroom(snakemake@input[[3]])
ref <- snakemake@config[["reference"]]
out <-  snakemake@output[[1]]

# mutations <- vroom("results/mutations/mutations_filtered.tsv")
# genes_gff <- readGff3("data/genomes/Fagus_sylvatica_v3.1.annot.gff", quiet = T)
# genes_tab <- vroom("data/genomes/Fagus_sylvatica_v3.1.annot.gff", skip = 1,
#                    col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))
# annotation <- vroom("data/genomes/Fagus_sylvatica_V3_genes_annotations_V3.csv")
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

overlaps <- interval_overlap(mutationsI, genes_gff)

genes <- lapply(1:length(overlaps), function(i){
  if(length(overlaps[[i]]) > 0)
    bind_cols(
      mutations[i,],
      genes_tab[overlaps[[i]],] %>% 
        dplyr::select(-seqid, -source) %>% 
        summarise(type = paste(type, collapse = ";"), start = min(start), end = max(end), score = dplyr::first(score),
                  strand = dplyr::first(strand), phase = dplyr::first(phase), attributes = paste(attributes, collapse = ";")))
  
}) %>% bind_rows() %>% 
  dplyr::select(SNV, type, start, end, score, strand, phase, attributes) %>% 
  dplyr::rename(gene_type = type, gene_start = start, gene_end = end, gene_score = score, gene_strand = strand,
                gene_phase = phase, gene_attributes = attributes)
  
genes <- genes %>% 
  separate(gene_attributes, "gene_id", ";") %>% 
  mutate(gene_id = gsub("ID=", "", gene_id)) %>% 
  left_join(annotation, by = c("gene_id"= "Fagus_assemblage_CEA_V3_geneID"))

write_tsv(genes, out)

