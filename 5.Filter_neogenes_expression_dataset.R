### 0. Libraries, functions and paths ------------------------------------------
library(readxl)         
library(dplyr)          
library(ggplot2)         
library(openxlsx)        
library(rtracklayer)     
library(GenomicRanges)   
library(Rgff)            
library(circlize)        
library(ComplexHeatmap)  

resdir <- ".../Neogenes/"
data_dir <- ".../2024-06-27/"
gtf_dir <- ".../2024-06-27/query_GTF_2/"

`%nin%` <- Negate(`%in%`)

### 1. Load data ---------------------------------------------------------------
Neogenes <- read_xlsx(paste0(data_dir, "EwNG_Neogene_2/means/EwNG_Neogene_2.all.means.xlsx"))
Neoantisens <- read_xlsx(paste0(data_dir, "EwNG_Neoantisens_2/means/EwNG_Neoantisens_2.all.means.xlsx"))
Putative <- read_xlsx(paste0(data_dir, "EwNG_Putative_Neogene_2/means/EwNG_Putative_Neogene_2.all.means.xlsx"))

final_annot <- read_xlsx(".../final.annot.unique.xlsx")
final_neo = final_annot %>%
  dplyr::filter(grepl("Ew_NG", Gene_curated))

final_neo = unique(final_neo$Gene_curated)

print(final_neo)

all <- rbind(Neogenes, Neoantisens)
all <- rbind(all, Putative)

names <- all$`0`
all = all[,2:ncol(all)]
rownames(all) <- names

### 2. Remove isoforms with too small sequences --------------------------------
# Function to calculate exon lengths for all genes
calculate_exon_lengths <- function(gtf_file) {
  # Import the GTF file
  gtf <- import(gtf_file)
  
  # Filter for exons
  exons <- gtf[gtf$type == "exon"]
  
  # Extract transcript_id from the attribute column
  transcript_ids <- sapply(strsplit(mcols(exons)$transcript_id, ";"), function(x) gsub('transcript_id "', '', x[1]))
  mcols(exons)$transcript_id <- transcript_ids
  
  # Calculate lengths of exons
  exon_lengths <- width(exons)
  
  # Create a data frame with transcript_id and exon_length
  exon_data <- data.frame(transcript_id = mcols(exons)$transcript_id, exon_length = exon_lengths)
  
  # Summarize lengths by transcript_id
  gene_lengths <- aggregate(exon_length ~ transcript_id, data = exon_data, sum)
  
  return(gene_lengths)
}


# Neogenes
gtf_file <- paste0(gtf_dir, "EwNG_Neogene.gtf")
Neogene_size <- calculate_exon_lengths(gtf_file)
print(Neogene_size)

# Neoantisens
gtf_file <- paste0(gtf_dir, "EwNG_Neoantisens.gtf")
Neoantisens_size <- calculate_exon_lengths(gtf_file)
print(Neoantisens_size)

# Other isoforms
gtf_file <- paste0(gtf_dir, "EwNG_Putative_Neogene.gtf")
Putative_size <- calculate_exon_lengths(gtf_file)
print(Putative_size)

all.length <- rbind(Neogene_size, Neoantisens_size)
all.length <- rbind(all.length, Putative_size)
all.length$exon_length = as.numeric(all.length$exon_length)

# Correlate mean expression in non-ewing datasets and query length fo each isoform 
all.2 <- as.data.frame(all)
all.2 = sapply(all.2, as.numeric)
all.2 = as.matrix(all.2)
rownames(all.2) = names
row_means <- rowMeans(all.2[, colnames(all.2) %nin% c("CUS-EW_illumina", "CUS-EW_nanopore")], na.rm = TRUE)
all.2 <- cbind(row_means , all.2)
FC_Ewing_Illumina <- all.2[,"CUS-EW_illumina"]/(all.2[,"row_means"]+0.0001)
FC_Ewing_Nanopore <- all.2[,"CUS-EW_nanopore"]/(all.2[,"row_means"]+0.0001)
all.2 <- cbind(FC_Ewing_Illumina , all.2)
all.2 <- cbind(FC_Ewing_Nanopore , all.2)
all.2 = as.data.frame(all.2)
all.2$transcript_id = names
all.2 = all.2 %>%
  dplyr::select(c(row_means, transcript_id, FC_Ewing_Illumina, FC_Ewing_Nanopore)) %>%
  arrange(transcript_id)
dim(all.2)
all.length = all.length %>%
  arrange(transcript_id)
dim(all.length)

size <- cbind(all.length, all.2)
size = size %>%
  dplyr::select(c(exon_length ,row_means, FC_Ewing_Illumina, FC_Ewing_Nanopore))

# Create a correlation plot between isoform expression and isoform length
plot <- ggplot(size, aes(x = log10(exon_length), y = row_means)) +
  geom_point(size = 3, alpha = 0.4, color = "darkblue") +  # Larger points with transparency and blue color
  geom_smooth(method = "lm", col = "darkblue") +  # Add a linear regression line 
  geom_vline(xintercept = 1.4, linetype = "dashed", color = "red", size = 1) +  # Add a vertical line at x=1.4 25bp
  theme_classic()+
  annotate("text", x = 1.4, y = 1000, label = "x=1.4 ~ 25bp", color = "black", size = 6, hjust = 0) +  # Annotation
  labs(title = "Correlation Plot",
       x = "log10(Exon Length)",
       y = "Row Means") +
  theme(axis.title.x = element_text(margin = margin(t = 10)),  # Increase distance between x-axis label and axis
        axis.title.y = element_text(margin = margin(r = 10)))  # Increase distance between y-axis label and axis

plot 

ggsave(plot=plot, filename = paste0(resdir, "correlation_expression_in_non_ewing_exon_length_neogenes.png"), units = "cm", width = 10, height = 10, dpi = 300)

### 3. Calculate fold change per gene Ewing vs non Ewing  ----------------------

# Calculate fold change per gene
original_rownames <- rownames(size)
new_rownames <- sub("\\..*", "", original_rownames)
size$gene<- new_rownames

size  <- size %>%
  group_by(gene) %>%
  summarize(across(everything(), mean, na.rm = TRUE)) %>%
  arrange(-FC_Ewing_Nanopore) %>%
  arrange(-FC_Ewing_Illumina)

rownames(size) <- size$gene

ord <- size$gene

# Neogenes not present in the expression dataframe because fully included in a known sequence
to_remove <- setdiff(final_neo, ord)
print(to_remove)

# Convert to matrix
names <- rownames(all)
all = sapply(all, as.numeric)
all = as.matrix(all)
rownames(all) = names

### 4. Prepare matrix expression per gene --------------------------------------
# Create gene annotation starting with isoform name
all.gene <- as.data.frame(all)
original_rownames <- rownames(all.gene)
new_rownames <- sub("\\..*", "", original_rownames)
all.gene$gene<- new_rownames

# Calculate mean expression per gene
gene_means <- all.gene %>%
  group_by(gene) %>%
  summarize(across(everything(), mean, na.rm = TRUE))

rownames(gene_means) = gene_means$gene
names <- rownames(gene_means)

# reorder matrix
ord = ord[ord %in% rownames(gene_means)]
gene_means <- gene_means[match(ord, rownames(gene_means)), ]
rownames(gene_means) = gene_means$gene
names = rownames(gene_means)
gene_means = gene_means %>% dplyr::select(-gene)
rownames(gene_means) = names

### 5. Remove genes expressed --------------------------------------------------
genes_expressed <- size %>% dplyr::filter(row_means>5)
dim(genes_expressed)
genes_expressed = rownames(genes_expressed)
print(genes_expressed)
dim(gene_means)
gene_means = gene_means %>% dplyr::filter(rownames(gene_means) %nin% genes_expressed)
dim(gene_means)

# Keep genes with high FC versus Ewing 
High_FC <- size %>% dplyr::filter(FC_Ewing_Illumina>15 | FC_Ewing_Nanopore>15)
High_FC = rownames(High_FC)
print(High_FC)
gene_means = gene_means %>% dplyr::filter(rownames(gene_means) %in% High_FC)
dim(gene_means)

# Retransform into matrix 
names = rownames(gene_means)
gene_means = sapply(gene_means, as.numeric)
gene_means = as.matrix(gene_means)
rownames(gene_means) = names

# Number of neogenes kept
dim(gene_means)

### 6. Plot heatmap ------------------------------------------------------------
# Color set up 
cols <-c("#6947A4", "#FBF6F4", "red")
mat = t(apply(gene_means, 1, scale))
colnames(mat) <- colnames(gene_means)

mat = t(mat)

# Extract and sort the row names
sorted_row_names <- sort(rownames(mat))

# Reorder the matrix based on sorted row names
ordered_mat <- mat[sorted_row_names, ]

colnames(ordered_mat)

# Neogenes excluded
expressed_neogenes <- setdiff(size$gene, colnames(ordered_mat))
print(expressed_neogenes)

# Add foldchange annotation 
fold_change_Illumina <- size
fold_change_Illumina = fold_change_Illumina[fold_change_Illumina$gene %in% rownames(gene_means),]
fold_change_Illumina = fold_change_Illumina$FC_Ewing_Illumina
names(fold_change_Illumina) <- rownames(gene_means)

fold_change_Nanopore <- size
fold_change_Nanopore = fold_change_Nanopore[fold_change_Nanopore$gene %in% rownames(gene_means),]
fold_change_Nanopore = fold_change_Nanopore$FC_Ewing_Nanopore
names(fold_change_Nanopore) <- rownames(gene_means)


# Create the annotation for fold change
fold_change_annotation <- HeatmapAnnotation(
  fold_change_Illumina = log10(0.0001 + fold_change_Illumina), 
  fold_change_Nanopore = log10(0.0001 + fold_change_Nanopore),
  col = list(
    fold_change_Illumina = colorRamp2(
      c(min(log10(0.0001 + fold_change_Illumina)), 0, max(log10(0.0001 + fold_change_Illumina))), 
      c("darkblue", "white", "red")
    ),
    fold_change_Nanopore = colorRamp2(
      c(min(log10(0.0001 + fold_change_Nanopore)), 0, max(log10(0.0001 + fold_change_Nanopore))), 
      c("darkblue", "white", "red")
    )
  ),
  annotation_height = unit(c(0.05, 0.05), "cm"),  # Adjust the height of the annotation bars
  annotation_legend_param = list(
    fold_change_Illumina = list(
      title = "log10(Fold Change Illumina)",
      title_gp = gpar(fontsize = 10, fontface="bold"),  # Title font size
      labels_gp = gpar(fontsize = 10), # Labels font size
      legend_height = unit(2.5, "cm"),   # Legend height
      legend_width = unit(1, "cm")     # Legend width
    ),
    fold_change_Nanopore = list(
      title = "log10(Fold Change Nanopore)",
      title_gp = gpar(fontsize = 10, fontface="bold"),  # Title font size
      labels_gp = gpar(fontsize = 10), # Labels font size
      legend_height = unit(2.5, "cm"),   # Legend height
      legend_width = unit(1, "cm")     # Legend width
    )
  )
)

cols <-c("#6947A4", "#FBF6F4", "red")

# Plot heatmap 
h <- Heatmap(
  ordered_mat, 
  col = cols, 
  name = "Expression Z-score",
  border_gp = gpar(col = "black", lty = 3),
  cluster_columns = FALSE,
  show_column_names = TRUE,
  cluster_rows = FALSE,
  show_row_names = TRUE,
  show_heatmap_legend = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  top_annotation = fold_change_annotation
)

# Draw the heatmap with the legend on the right side
draw(h, heatmap_legend_side = "right")

png(paste0(resdir, "heatmap_neogenes_per_gene_FC15_Rowmeans_5.png"), width = 28, height = 18, units = "cm", res = 300)
draw(h)
dev.off()

### 7. Create new gtf and dataframe without neogenes expessed in public datasets ---
  # 7.1 Export excel table -----------------------------------------------------

neogenes_to_remove <- c(to_remove, expressed_neogenes)
print(neogenes_to_remove)

dim(final_annot)
final_annot = final_annot %>%
  dplyr::filter(Gene_curated %nin% neogenes_to_remove) %>%
  dplyr::filter(grepl("Ew_NG", Gene_curated)) %>%
  dplyr::mutate(Gene_curated = ifelse(Gene_curated == "Putative_Ew_NG69", "Ew_NG69", Gene_curated),
                isoform_type = ifelse(Gene_curated == "Putative_Ew_NG69", "Neogene", isoform_type),
                isoform_version = ifelse(isoform_version == "Putative_Ew_NG69.1", "Ew_NG69.1", isoform_version))
  

dim(final_annot)
  
write.xlsx(final_annot, file = ".../Neogenes_10_07_2024.xlsx", sep="\t")

  # 7.2 Export gtf -------------------------------------------------------------
# Load gtf 
gtf <- readGFF(".../Novel_isoforms_final.gtf")
dim(gtf)
class(gtf)

# Remove columns not used 
gtf.novel = gtf %>%
  dplyr::filter(gene_id %nin% neogenes_to_remove) %>%
  dplyr::filter(isoform_type %in% c("Neogene", "Neoantisens", "Putative_Neogene")) %>%  # Keep only neogenes
  dplyr::mutate(gene_id = ifelse(gene_id == "Putative_Ew_NG69", "Ew_NG69", gene_id),
                gene_name = ifelse(gene_name == "Putative_Ew_NG69", "Ew_NG69", gene_name),
                isoform_type = ifelse(gene_id == "Ew_NG69", "Neogene", isoform_type),
                transcript_id = ifelse(transcript_id == "Putative_Ew_NG69.1", "Ew_NG69.1", transcript_id))

dim(gtf.novel)

unique(gtf.novel$gene_id)
unique(gtf.novel$gene_name)
unique(gtf.novel$transcript_id)

# Remove neogenes expressed in public datasets ----------------------------

gtf.novel = makeGRangesFromDataFrame(gtf.novel, keep.extra.columns=TRUE,
                                     seqinfo=NULL)

rtracklayer::export(gtf.novel, ".../Neogenes_10072024.gtf")


