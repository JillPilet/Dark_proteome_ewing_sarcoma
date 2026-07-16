# Load libraries----------------------------------------------------------------

library(readxl)
library(dplyr)    
library(openxlsx)  
`%nin%` <- Negate(`%in%`)

# Load annotation  and modify protein annotation -------------------------------
Neoisoforms <- readxl::read_xlsx(".../Neoisoforms_10_07_2024.xlsx")
Neogenes <- readxl::read_xlsx(".../Neogenes_10_07_2024.xlsx")

annot <- rbind(Neoisoforms, Neogenes)

# Re-annotate isoform type using the truncated isoforms manual curation --------
manual_curation <- readxl::read_xlsx(".../Truncated_isoforms_manual_annot.xlsx")

curated_truncated = manual_curation %>% filter(Manual_curation_truncated=="Truncated_isoform")
curated_truncated = curated_truncated$Gene_automatic_annot

# Now create a new annotation of isoforms---------------------------------------
annot = annot %>%
  mutate(isoform_recode = case_when(isoform_type =="Neogene" ~ "Neogene", 
                                    isoform_type =="Neoantisens" ~ "Neogene", 
                                    isoform_type =="Long_isoform" ~ "Non-canonical", 
                                    isoform_type =="Neoisoform" ~ "Non-canonical",
                                    isoform_type =="Truncated_isoform" ~ "Truncated_isoform"),
         isoform_final = case_when(Gene_curated %in%  curated_truncated ~ isoform_recode,
                                   isoform_recode =="Truncated_isoform" ~ "Non-canonical", 
                                   isoform_recode =="Non-canonical" ~ "Non-canonical", 
                                   isoform_recode =="Neogene" ~ "Neogene"))

table(annot$isoform_final)

# Annotate isoforms with neoexons ----------------------------------------------
Neoexons <- readxl::read_xlsx(".../Neoexons_10_07_2024.xlsx")
dim(Neoexons)
Neoexons = Neoexons %>% dplyr::select(c("isoform_version", "Gene_curated"))
dim(Neoexons)
Neoexons$Neoexon = "Yes"

annot = left_join(annot, Neoexons, by="isoform_version")

table(annot$isoform_final, annot$Neoexon)

annot = annot %>%
  mutate(isoform_final = case_when(Neoexon =="Yes" & isoform_final=="Truncated_isoform" ~ "Neoexon; Intragenic TSS isoform", 
                                   Neoexon =="Yes" ~ "Neoexon",
                                   isoform_final=="Truncated_isoform" ~ "Intragenic TSS isoform", 
                                   isoform_final =="Non-canonical" ~ "Non-canonical", 
                                   isoform_final =="Neogene" ~ "Neogene"))

table(annot$isoform_final)

# Export isoform complete annotation -------------------------------------------
write.xlsx(annot, ".../neogenes_isoforms_annot_20260320.xlsx")
