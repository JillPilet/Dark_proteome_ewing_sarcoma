# Load libraries ---------------------------------------------------------------
library(readxl)
library(dplyr)
library(openxlsx)

# Read table -------------------------------------------------------------------
curie_neo <- readxl::read_xlsx(".../neoisoforms_ORFs_with_dups_collasped.xlsx")

# Look at the table ------------------------------------------------------------
dim(curie_neo)

length(unique(curie_neo$isoform))
length(unique(curie_neo$Gene_curated))

table(curie_neo$isoform_final)

# Select columns of interest ---------------------------------------------------

curie_neo_clean = curie_neo %>%
  dplyr::select(c(Gene_curated, isoform_version, merged_ids,
                  TSS, Strand, # Isoforms coordinates
                  Detected.in.A673, Detected.in.TC71, Detected.in.EW7, Detected.in.EW16, Detected.in.Ewima1, # PacBio & Nanopore detection at isoform level
                  Detected.in.Nanopore, Detected.in.PacBio,
                  isoform_final, isoform_final_subtype, structural_category.x, subcategory, # isoform sqanti annotation
                  FLI1_GGAAm_in_TAD, TAD_name, A673_log2FC_High_Low, TC71_72h_log2FC_High_Low, Median_Transcript_MSC.6., Nb_reads_Nanopore, # Filters to get EwS specific neoisoforms and neogenes 
                  neotranscript_orf_type_per_isoform, known_orf_id, ORF_ranges, Protein
                  ))


# Rename columns of interest ---------------------------------------------------

curie_neo_clean = curie_neo_clean %>%
  mutate(Detected.in.A673 = ifelse(Detected.in.A673==0, "No", "Yes"),
         Detected.in.TC71 = ifelse(Detected.in.TC71==0, "No", "Yes"),
         Detected.in.EW7 = ifelse(Detected.in.EW7==0, "No", "Yes"),
         Detected.in.EW16 = ifelse(Detected.in.EW16==0, "No", "Yes"),
         Detected.in.Ewima1 = ifelse(Detected.in.Ewima1==0, "No", "Yes"),
         Detected.in.Nanopore = ifelse(Detected.in.Nanopore==0, "No", "Yes"),
         Detected.in.PacBio = ifelse(Detected.in.PacBio==0, "No", "Yes")) %>%                        
  dplyr::rename(Detected_in_A673 = Detected.in.A673,
                Detected_in_TC71 = Detected.in.TC71,
                Detected_in_EW7 = Detected.in.EW7,
                Detected_in_EW16 = Detected.in.EW16,
                Detected_in_Ewima1 = Detected.in.Ewima1,
                Detected_in_Nanopore = Detected.in.Nanopore,
                Detected_in_PacBio = Detected.in.PacBio,
                Neoisoform = isoform_version,
                `Neoisoform_ID_in_each_cell_line`= merged_ids,
                `Gene_name` = Gene_curated,
                `Neoisoform_annotation` = isoform_final,
                `Neoisoform_annotation_subtype` = isoform_final_subtype,
                `Sqanti_structural_category` = structural_category.x,
                `Sqanti_structural_subcategory`= subcategory,
                `Median_Transcript_expression_in_6_MSCs` = Median_Transcript_MSC.6.,
                `Number_of_Nanopore_reads` = Nb_reads_Nanopore,
                `Neoisoform_ORF_category` = neotranscript_orf_type_per_isoform,
                `ORFs_ids`= known_orf_id,
                `ORF_ranges` = ORF_ranges,
                `Protein_sequence` = Protein)
  
  
write.xlsx(curie_neo_clean, ".../Supplementary_table_S6.xlsx")
