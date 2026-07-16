### 0. Libraries, functions and paths ------------------------------------------
library(readxl)         
library(dplyr)           
library(ComplexHeatmap) 
library(extrafont)       
loadfonts()

resdir <- ".../Novel_exons/"
data_dir <- ".../Expression_PM/"

`%nin%` <- Negate(`%in%`)

### 1. Load data ---------------------------------------------------------------
Meta_data <- read.csv(paste0(data_dir, "EwS_ribo_meta_all.csv")) %>%
  dplyr::filter(source == "cell_line", type != "A673_scrambl") %>%
  arrange(type, condition_dox)

all <- read.csv(paste0(data_dir, "ppm_counts_all.csv"))
ORF_ids <- read.csv(".../neoisoform_orf_table.csv")

ORF_annot_isoforms <- readxl::read_xlsx(".../neoisoforms_ORFs_with_dups.xlsx")

ORF_annot_isoforms = ORF_annot_isoforms %>% 
  as.data.frame() %>%
  dplyr::select(known_orf_id, isoform_final) %>%
  dplyr::filter(known_orf_id!="NA")

# Modify annotations 
ORF_ids = left_join(ORF_ids, ORF_annot_isoforms, by="known_orf_id")

table(ORF_ids$neotranscript_orf_type.x)
table(ORF_ids$isoform_final)

ORF_ids = ORF_ids %>%
  mutate(neotranscript_orf_type_final = case_when(neotranscript_orf_type.x=="In-frame non-canonical" ~ "In-frame ORF",
                                                  neotranscript_orf_type.x=="In-frame truncated" ~ "In-frame ORF",
                                                  neotranscript_orf_type.x=="Out-of-frame novel" ~ "Out-of-frame ORF"))
table(ORF_ids$neotranscript_orf_type_final)
ORF_ids_to_keep <- ORF_ids$known_orf_id

# Keep only isoforms to plot
all <- all %>% dplyr::filter(rownames(.) %in% ORF_ids_to_keep)

### 2. Prepare metadata and expression matrix ----------------------------------
Meta_data <- Meta_data %>%
  mutate(id = case_when(
    source == "cell_line" ~ paste0("CL_", sample_id),
    source == "patient_sample" ~ paste0("PS_", sample_id)
  ))

# Keep only samples present in metadata
names <- rownames(all)
all <- all[, colnames(all) %in% Meta_data$id]
all <- as.matrix(sapply(all, as.numeric))
rownames(all) <- names

# Align columns with metadata
all <- all[, match(Meta_data$id, colnames(all))]
stopifnot(all(colnames(all) == Meta_data$id))

### 3. Normalize each cell line relative to its -DOX condition -----------------
all_norm <- matrix(NA, nrow=nrow(all), ncol=ncol(all))
rownames(all_norm) <- names
colnames(all_norm) <- colnames(all)

cell_lines <- unique(Meta_data$type)

for(cl in cell_lines){
  samples_cl <- Meta_data$id[Meta_data$type == cl]
  samples_neg <- Meta_data$id[Meta_data$type == cl & Meta_data$condition_dox == "neg"]
  
  # Mean expression of -DOX isoforms
  baseline <- rowMeans(all[, samples_neg, drop=FALSE], na.rm=TRUE)
  
  # Subtract baseline
  all_norm[, samples_cl] <- sweep(all[, samples_cl, drop=FALSE], 1, baseline, FUN="-")
}

# Optional: Z-score after -DOX normalization
mat <- t(apply(all_norm, 1, scale))
colnames(mat) <- colnames(all_norm)
mat <- t(mat)

### 4. Merge ORF category for row ordering -------------------------------------
ORF_ids_sub <- ORF_ids[ORF_ids$known_orf_id %in% colnames(mat),
                       c("known_orf_id","neotranscript_orf_type_final", "isoform_final")]

ORF_ids_sub <- ORF_ids_sub[match(colnames(mat), ORF_ids_sub$known_orf_id), ]

# Define category order
category_order <- c("In-frame ORF", "Out-of-frame ORF")
ORF_ids_sub$neotranscript_orf_type_final <- factor(ORF_ids_sub$neotranscript_orf_type_final,
                                               levels=category_order)

# Sort rows by category
sorted_idx <- order(ORF_ids_sub$neotranscript_orf_type_final, ORF_ids_sub$isoform_final)
mat_sorted <- mat[,sorted_idx]
ORF_ids_sorted <- ORF_ids_sub[sorted_idx, ]

### 5. Column annotation -------------------------------------------------------
Meta_data$condition_dox[Meta_data$source != "cell_line"] <- NA

type <- c("A673"="#fdcf9d", "MHHES"="#00AAFF", "RDES" = "lightblue", "SKNMC" = "green", "TC106" = "yellow")
dox_col <- c("pos"="#b4b4fa", "neg"="gray90")


ha <- HeatmapAnnotation(
  type = Meta_data$type,
  Dox = Meta_data$condition_dox,
  col = list(type = type, Dox = dox_col),
  na_col="white",
  annotation_name_gp = gpar(fontsize=10, fontfamily="Helvetica"),
  simple_anno_size = unit(5, "mm"),
  annotation_legend_param = list(
    type=list(title="Sample type", title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=9)),
    Dox=list(title="Doxycycline", title_gp=gpar(fontsize=10), labels_gp=gpar(fontsize=9))
  )
)

### 6. Row annotation (category) ------------------------------------------------
row_ha <- rowAnnotation(
  Category = ORF_ids_sorted$neotranscript_orf_type_final,
  isoform = ORF_ids_sorted$isoform_final,
  col = list(Category = c("In-frame ORF"="#6EA9D2",
                          "Out-of-frame ORF"="#D0E8F8"),
             isoform = c("Neoexon" = "#EA7580FF", 
                         "Neoexon; Intragenic TSS isoform" = "#088BBEFF",
                         "Intragenic TSS isoform"= "#F8CD9CFF",
                         "Non-canonical"= "#1BB6AFFF")),
  show_annotation_name=TRUE
)

### 7. Heatmap colors ----------------------------------------------------------
cols <- c("blue", "white", "red")

### 8. Create heatmap ----------------------------------------------------------
flipped_mat <- t(mat_sorted)

h <- Heatmap(
  flipped_mat,
  col = cols,
  name = "Expression Z-score (-DOX normalized)",
  border_gp=gpar(col="black", lty=3),
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  show_row_names=FALSE,
  show_row_dend = FALSE,
  show_column_names=T,
  column_names_gp=gpar(fontsize=12, fontfamily="Helvetica"),
  row_names_gp=gpar(fontsize=12, fontfamily="Helvetica"),
  heatmap_legend_param=list(
    title_gp=gpar(fontsize=12, fontfamily="Helvetica", fontface="plain"),
    labels_gp=gpar(fontsize=12, fontfamily="Helvetica"),
    legend_height=unit(3,"cm"),
    legend_width=unit(4,"cm")
  ),
  top_annotation=ha,
  left_annotation=row_ha
)

draw(h, heatmap_legend_side="right", annotation_legend_side="right")

### 9. Save to PDF -------------------------------------------------------------
pdf(paste0(resdir,"heatmap_all_isoforms_DOX_norm_by_category.pdf"), width=6, height=4, fonts="Helvetica")
draw(h)
dev.off()

