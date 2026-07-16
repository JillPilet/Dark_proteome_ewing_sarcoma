# Load packages
library(rtracklayer)
library(openxlsx)
library(dplyr)
library(readxl)

# Load data
final_annot <-  read_xlsx("C:/Users/Administrateur/Documents/Neogene_long_read_project/All_merged/V9/Merged_file/Neoisoforms_10_07_2024.xlsx", sep="\t")
gtf <- readGFF("C:/Users/Administrateur/Documents/Neogene_long_read_project/All_merged/V9/gtf/Novel_isoforms_final.gtf")
dim(gtf)
class(gtf)

gtf.novel = gtf %>%
  dplyr::filter(transcript_id %in% final_annot$isoform_version)
dim(gtf.novel)

unique(gtf.novel$gene_id)
unique(gtf.novel$transcript_id)
length(unique(gtf.novel$gene_id))
length(unique(gtf.novel$transcript_id))

gtf.novel = makeGRangesFromDataFrame(gtf.novel, keep.extra.columns=TRUE,
                                     seqinfo=NULL)

rtracklayer::export(gtf.novel, "C:/Users/Administrateur/Documents/Neogene_long_read_project/All_merged/V9/gtf/Neoisoforms_10_07_2024.gtf")
