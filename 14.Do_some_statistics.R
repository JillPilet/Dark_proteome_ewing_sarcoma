# The goal of this script is to count isoforms

# Load libraries ---------------------------------------------------------------
library(readxl) 
library(dplyr)    

# Load tables ------------------------------------------------------------------

neo_with_dup <- readxl::read_xlsx(".../neoisoforms_ORFs_with_dups.xlsx")

nb_isoform_with_translon_dup = neo_with_dup %>% filter(neotranscript_orf_type!="No ORF") 
dim(nb_isoform_with_translon_dup)
nb_isoform_with_translon_dup$neotranscript_orf_type
length(unique(nb_isoform_with_translon_dup$isoform_version))  

unique(nb_isoform_with_translon_dup$Gene_curated)

# Neoexons ---------------------------------------------------------------------
table(neo_with_dup$isoform_final)
neoexons = neo_with_dup %>% filter(isoform_final %in% c("Neoexon", "Neoexon; Intragenic TSS isoform"))
length(unique(neoexons$isoform_version))
neoexon_translon =  neoexons %>% filter(neotranscript_orf_type!="No ORF") 
length(unique(neoexon_translon$isoform_version))
unique(neoexon_translon$Gene_curated)

neoexon_translon_unique = neoexon_translon %>% filter(!duplicated(known_orf_id))
table(neoexon_translon_unique$neotranscript_orf_type)

# Truncated---------------------------------------------------------------------
Truncated = neo_with_dup %>% filter(isoform_final %in% c("Intragenic TSS isoform", "Neoexon; Intragenic TSS isoform"))
length(unique(Truncated$isoform_version))
length(unique(Truncated$Gene_curated))

# non canonical isoforms--------------------------------------------------------
nc = neo_with_dup %>% filter(isoform_final %in% c("Non-canonical"))
length(unique(nc$isoform_version))
nc_translons = nc %>% filter(neotranscript_orf_type!="No ORF") 
nc_translons_unique = nc_translons %>% filter(!duplicated(known_orf_id))
table(nc_translons_unique$neotranscript_orf_type)

# Neoexon + truncated ----------------------------------------------------------
neoexon_truncated = neo_with_dup %>% filter(isoform_final %in% c("Neoexon", 
                                                                 "Intragenic TSS isoform",
                                                                 "Neoexon; Intragenic TSS isoform"))

neoexon_truncated_translon =  neoexon_truncated %>% filter(neotranscript_orf_type!="No ORF") 
neoexon_truncated_translon_unique = neoexon_truncated_translon %>% filter(!duplicated(known_orf_id))
table(neoexon_truncated_translon_unique$neotranscript_orf_type)

