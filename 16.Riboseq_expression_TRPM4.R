# Libraries, functions and paths -----------------------------------------------
library(dplyr)
resdir <- ".../Novel_exons/"
data_dir <- ".../Expression_PM/"

`%nin%` <- Negate(`%in%`)

# Load data --------------------------------------------------------------------
Meta_data <- read.csv(paste0(data_dir, "EwS_ribo_meta_all.csv")) %>%
  dplyr::filter(source == "patient_sample") 

all <- read.csv(paste0(data_dir, "ppm_counts_all.csv"))

# Keep only TRPM4 ORFs and dORF-------------------------------------------------
canonical_TRPM4 <- rownames(all)[grep("ENSG00000130529", rownames(all))]
canonical_TRPM4

all <- all %>% dplyr::filter(rownames(.) %in% c(canonical_TRPM4, "TCONS_00433586_5495_5734"))

# Keep only samples present in metadata-----------------------------------------
Meta_data <- Meta_data %>%
  mutate(id = case_when(
    source == "cell_line" ~ paste0("CL_", sample_id),
    source == "patient_sample" ~ paste0("PS_", sample_id)
  ))


names <- rownames(all)
all <- all[, colnames(all) %in% Meta_data$id]
all <- as.matrix(sapply(all, as.numeric))
rownames(all) <- names

all <- all %>%
  as.data.frame() %>%
  mutate(
    type = ifelse(rownames(.) == "TCONS_00433586_5495_5734", "dORF", "canonical"),
    row_mean = rowMeans(across(where(is.numeric)), na.rm = TRUE)
  )

dORF_mean <- all %>%
  dplyr::filter(type == "dORF") %>%
  dplyr::pull(row_mean)
dORF_mean

canonical_mean <- all %>%
  dplyr::filter(type == "canonical") %>%
  summarise(mean(row_mean, na.rm = TRUE)) %>%
  pull()

fold_change <- dORF_mean / canonical_mean
fold_change
