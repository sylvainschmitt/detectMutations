library(tidyverse)
library(vroom)

mutations <- vroom(snakemake@input[["mutations"]])
validation <- vroom(snakemake@input[["validation"]])
te <- vroom(snakemake@input[["te"]])
genes <- vroom(snakemake@input[["genes"]])
spectra <- vroom(snakemake@input[["spectra"]])
ref <- snakemake@config[["reference"]]
out <-  snakemake@output[[1]]

# mutations <- vroom("results/mutations/mutations_filtered.tsv")
# validation <- vroom("results/mutations/cross_validation.tsv")
# spectra <- vroom("results/mutations/spectra.tsv")
# genes <- vroom("results/mutations/genes.tsv")
# te <- vroom("results/mutations/te.tsv")
# ref <- "Fagus_sylvatica_v3"

mutations <- mutations %>% 
  filter(reference == ref) %>% 
  mutate(SNV = paste0(CHROM, "#", POS))

validation <- validation %>%
  reshape2::dcast(SNV ~ reference) %>% 
  dplyr::select(-Fagus_sylvatica_v3) %>% 
  mutate(validation  = "all") %>% 
  mutate(validation = ifelse(Fagus_sylvatica_mutant_v0 == 0, "only revertant", validation)) %>% 
  mutate(validation = ifelse(Fagus_sylvatica_revertant_v0 == 0, "only mutant", validation)) %>% 
  dplyr::select(SNV, validation)

mutations <- left_join(mutations, validation)

mutations <- left_join(mutations, spectra)

mutations <- left_join(mutations, genes)

mutations <- left_join(mutations, te) %>% 
  separate(te_attributes, "te_id", ";") %>% 
  mutate(te_id = gsub("ID=", "", te_id))

mutations <- mutations %>% 
  dplyr::rename(scaffold = CHROM, position = POS, ref = REF, alt = ALT, strelka2_filter = FILTER,
                qss = QSS, tqss = TQSS, qss_nt = QSS_NT, tqss_nt = TQSS_NT, dp = DP, mq = MQ, mq0 = MQ0, 
                read_position_rank_sum = ReadPosRankSum, snvsb = SNVSB, evs = SomaticEVS,
                mutation_dp = mutation_DP, mutation_fdp = mutation_FDP, mutation_sdp = mutation_SDP,
                mutation_subdp = mutation_SUBDP, mutation_alt_count_t1 = mutation_altCountT1, mutation_alt_count_t2 = mutation_altCountT2,
                mutation_ref_count_t1 = mutation_refCountT1, mutation_ref_count_t2 = mutation_refCountT2, mutation_af = mutation_AF,
                normal_dp = normal_DP, normal_fdp = normal_FDP, normal_sdp = normal_SDP,
                normal_subdp = normal_SUBDP, normal_alt_count_t1 = normal_altCountT1, normal_alt_count_t2 = normal_altCountT2,
                normal_ref_count_t1 = normal_refCountT1, normal_ref_count_t2 = normal_refCountT2, normal_af = normal_AF,
                tumor = tumor, normal = normal, reference = reference, mutation_filter = Filter, snv = SNV,
                cross_validation = validation, mutation_class = class, mutation_type = type, mutation_spectra = spectra,
                gene_type = gene_type, gene_start = gene_start, gene_end = gene_end, gene_score = gene_score, gene_strand = gene_strand,
                gene_phase = gene_phase, gene_id = gene_id, annotation_defline = Defline, gene_id_Mishra2018 = geneID_genome_hetre_Mishra2018,
                tair10_gene_id = TAIR10_geneID, tair10_go = TAIR10_GO, tair10_function = TAIR10_function, 
                quercus_haplov2.3_gene_id = Quercus_haploV2.3_geneID, quercus_haplov2.3_function = Quercus_haploV2.3_funcat,
                iprscan_go = iprscan_GO, kog = KOG, ec_id = `EC ID`, panther_id = PANTHER_ID, pfam_id = `PFAM ID`, 
                te_type = te_type, te_start = te_start, te_end = te_end, te_score = te_score, te_strand = te_strand,
                te_phase = te_phase, te_id = te_id) %>% 
  dplyr::select(
    snv, reference, tumor, normal, scaffold, position, ref, alt, # mutation characteristics
    cross_validation, mutation_filter, strelka2_filter, # mutations filters
    qss, tqss, qss_nt, tqss_nt, dp, mq, mq0, read_position_rank_sum, snvsb, evs, # mutations general properties
    mutation_dp, mutation_fdp, mutation_sdp, mutation_subdp, mutation_alt_count_t1, # mutated tissue properties
    mutation_alt_count_t2, mutation_ref_count_t1, mutation_ref_count_t2, mutation_af, # mutated tissue properties
    normal_dp, normal_fdp, normal_sdp, normal_subdp, normal_alt_count_t1, # normal tissue properties
    normal_alt_count_t2, normal_ref_count_t1, normal_ref_count_t2, normal_af, # normal tissue properties
    mutation_class, mutation_type, mutation_spectra, # mutations spectra
    gene_id, gene_type, gene_start, gene_end, gene_score, gene_strand, gene_phase, # genes
    annotation_defline, gene_id_Mishra2018, tair10_gene_id, tair10_go, tair10_function, # gene annotations
    quercus_haplov2.3_gene_id, quercus_haplov2.3_function, iprscan_go, kog, ec_id, panther_id, pfam_id, # gene annotations
    te_id, te_type, te_start, te_end, te_score, te_strand, te_phase
  )





vroom_write(mutations, out)