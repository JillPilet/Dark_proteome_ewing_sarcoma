# Load libraries----------------------------------------------------------------

library(readxl)
library(dplyr)
library(rtracklayer)
library(tidyr)
library(Rgff)
library(GenomicRanges)
library(openxlsx)

# Load annotation and gtf ------------------------------------------------------
annot <- read_xlsx(".../final.annot.xlsx")
gtf <- readGFF(".../subset_unique_neoisoforms.gtf")

`%nin%` <- Negate(`%in%`)

# Add isoform version ----------------------------------------------------------
annot <- annot %>%
  filter(!duplicated(isoform)) %>%
  group_by(Gene_curated) %>%
  mutate(row_num = row_number())

annot$isoform_version <- paste(annot$Gene_curated, annot$row_num, sep = ".")

# Select columns of interest 
short.annot <- annot %>%
  select(c(isoform, Gene_curated, isoform_version, isoform_type, structural_category, Gene_Annot_if_different, associated_transcript, membrane_annot, Gene_automatic_annot))

# Splitting Gene_Annot_if_different column into two new columns
short.annot <- separate(short.annot, Gene_Annot_if_different, into = c("gffcompare_Ensembl_gene_id", "gffcompare_Ensembl_transcript_id"), sep = "\\|")

colnames(short.annot) <- c("transcript_id", "Gene_curated", "isoform_version", "isoform_type", "structural_category", "gffcompare_Ensembl_gene_id", "gffcompare_Ensembl_transcript_id", "associated_transcript", "membrane_annot", "Gene_gff_compare_old_annot")

in_annot_only = short.annot[short.annot$transcript_id %nin% gtf$transcript_id,]
in_gtf_only = gtf[gtf$transcript_id %nin% short.annot$transcript_id ,]

dim(in_annot_only)
length(in_gtf_only$gene_name)

gtf.novel <- left_join(gtf, short.annot, by="transcript_id")

gtf.novel$merge_id = gtf.novel$transcript_id
gtf.novel$XLOC_id = gtf.novel$gene_id

dim(gtf.novel)

# Clean aannotations
gtf.novel = gtf.novel %>%
  select(!transcript_id) %>%
  select(!gene_id) %>%
  select(!gene_name) %>%
  dplyr::rename(., transcript_id = isoform_version,
         gene_id = Gene_curated) %>%
  mutate(gene_name = ifelse(type=="transcript", paste0(gene_id), paste0(NA))) %>%
  filter(merge_id %in% annot$isoform) %>%
  mutate(gffcompare_Ensembl_transcript_id = ifelse(isoform_type %in% c("Neogene", "Neoantisens", "Putative_Neogene"), "na", gffcompare_Ensembl_transcript_id)) %>%
  mutate(gffcompare_Ensembl_gene_id = ifelse(isoform_type %in% c("Neogene", "Neoantisens", "Putative_Neogene"), "na", gffcompare_Ensembl_gene_id)) %>%
  mutate(gff_compare_sqanti_Ensembl_transcript_id = ifelse(is.na(gffcompare_Ensembl_transcript_id), associated_transcript, gffcompare_Ensembl_transcript_id))

gtf.novel[gtf.novel$merge_id %nin% short.annot$transcript_id ,]

dim(gtf.novel)

# Select columns of interest for export
gtf.novel = gtf.novel[,c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "gene_id", "transcript_id", "gene_name", "oId",
             "cmp_ref","class_code","tss_id","num_samples","exon_number","contained_in","cmp_ref_gene", "merge_id", "XLOC_id","isoform_type", "structural_category", "gffcompare_Ensembl_gene_id", "gffcompare_Ensembl_transcript_id", "associated_transcript", "gff_compare_sqanti_Ensembl_transcript_id", "membrane_annot","Gene_gff_compare_old_annot")]


# Make GRanges object
gtf.novel = makeGRangesFromDataFrame(gtf.novel, keep.extra.columns=TRUE,
                                     seqinfo=NULL)

# Export gtf and corresponding table
rtracklayer::export(gtf.novel, ".../Novel_isoforms_final.gtf")

write.xlsx(annot, file = ".../final_annot_unique.xlsx", sep="\t")

