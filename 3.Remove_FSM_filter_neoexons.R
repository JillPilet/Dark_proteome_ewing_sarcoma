# The goal of this script is to remove full-splice match isoform and annotate isoforms with neoexon

### 0. Libraries, functions and paths ------------------------------------------
library(ggplot2)         
library(dplyr)           
library(readxl)          
library(openxlsx)        
library(rtracklayer)    
library(GenomicRanges)   
library(Rgff)            
library(circlize)        
library(ComplexHeatmap)  
library(paletteer)       
library(extrafont)  

resdir <- ".../Novel_exons/"
data_dir <- ".../Quantification_public_datasets_Erkan/2024-06-27/"
gtf_dir <- ".../Quantification_public_datasets_Erkan/2024-06-27/query_GTF_1/"

`%nin%` <- Negate(`%in%`)

### 1. Load data ---------------------------------------------------------------
final_annot <- read_xlsx(".../final.annot.unique.xlsx")
Truncated <- read_xlsx(paste0(data_dir, "EwNG_Truncated_isoform_1/means/EwNG_Truncated_isoform_1.all.means.xlsx"))
Long <- read_xlsx(paste0(data_dir, "EwNG_Long_isoform_1/means/EwNG_Long_isoform_1.all.means.xlsx"))
Others <- read_xlsx(paste0(data_dir, "EwNG_Neoisoform_1/means/EwNG_Neoisoform_1.all.means.xlsx"))

all <- rbind(Truncated, Long)
all <- rbind(all, Others)

names <- all$`0`
all = all[,2:ncol(all)]
rownames(all) <- names

### 2. Filter out known isoforms -----------------------------------------------
# Restrict annotation to novel spliced forms 
final_annot = final_annot %>%
  dplyr::filter(isoform_type %in% c("Truncated_isoform", "Neoisoform", "Long_isoform"))

dim(final_annot)
length(unique(final_annot$Gene_curated))

final_annot = final_annot %>%
  dplyr::filter(structural_category !="full-splice_match")
dim(final_annot)
length(unique(final_annot$Gene_curated))

table(final_annot$ENSG_class_code) # No "=" class

#write.xlsx(final_annot, file = ".../Neoisoforms_10_07_2024.xlsx", sep="\t")

# Extract the relevant part of the rownames
novel_rownames <- gsub("_.*", "", rownames(all))
rownames(all) = novel_rownames
all = all %>% dplyr::filter(rownames(.) %in% final_annot$isoform_version)
dim(all)

count_genes = final_annot %>% dplyr::filter(isoform_version %in% rownames(all))
length(unique(count_genes$isoform_version))
length(unique(count_genes$Gene_curated))

### 3. Remove small queries ----------------------------------------------------
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


# Long isoforms
gtf_file <- paste0(gtf_dir, "EwNG_Long_isoform.gtf")
Long_isoforms_size <- calculate_exon_lengths(gtf_file)
print(Long_isoforms_size)

# Truncated isoforms
gtf_file <- paste0(gtf_dir, "EwNG_Truncated_isoform.gtf")
Truncated_isoforms_size <- calculate_exon_lengths(gtf_file)
print(Truncated_isoforms_size)

# Other isoforms
gtf_file <- paste0(gtf_dir, "EwNG_Neoisoform.gtf")
Other_isoforms_size <- calculate_exon_lengths(gtf_file)
print(Other_isoforms_size)

all.length <- rbind(Long_isoforms_size, Truncated_isoforms_size)
all.length <- rbind(all.length, Other_isoforms_size)
all.length$exon_length = as.numeric(all.length$exon_length)

novel_transcript_id <- gsub("_.*", "", all.length$transcript_id)
all.length$transcript_id = novel_transcript_id
all.length = all.length %>% dplyr::filter(transcript_id %in% rownames(all))
dim(all.length)

# Correlate mean expression in non-ewing datasets and query length 
names= rownames(all)
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

# Create a correlation plot
plot <- ggplot(size, aes(x = log10(exon_length), y = row_means)) +
  geom_point(size = 3, alpha = 0.4, color = "darkblue") +  # Larger points with transparency and blue color
  geom_smooth(method = "lm", col = "darkblue") +  # Add a linear regression line 
  geom_vline(xintercept = 1.4, linetype = "dashed", color = "red", size = 1) +  # Add a vertical line at x=1.4 25bp
  theme_classic()+
  annotate("text", x = 1.4, y = 20000, label = "x=1.4 ~ 25bp", color = "black", size = 6, hjust = 0) +  # Annotation
  labs(title = "Correlation Plot",
       x = "log10(Exon Length)",
       y = "Row Means") +
  theme(axis.title.x = element_text(margin = margin(t = 10)),  # Increase distance between x-axis label and axis
        axis.title.y = element_text(margin = margin(r = 10)))  # Increase distance between y-axis label and axis

plot 

ggsave(plot=plot, filename = paste0(resdir, "correlation_expression_in_non_ewing_exon_length.png"), units = "cm", width = 10, height = 10, dpi = 300)

# reorder matrix
size  <- size %>%
  arrange(-FC_Ewing_Nanopore) %>%
  arrange(-FC_Ewing_Illumina)

ord <- rownames(size)

ord = ord[ord %in% rownames(all)]
all <- all[match(ord, rownames(all)), ]
rownames(all) = ord
names = rownames(all)

# Remove small queries
small_query <- size %>% dplyr::filter(exon_length<=25)
dim(small_query)
small_query = rownames(small_query)
dim(all)
all = all %>% dplyr::filter(rownames(all) %nin% small_query)
dim(all)

nb_genes <-  gsub("\\.[0-9]+$", "", rownames(all))
unique(nb_genes)
length(unique(nb_genes))


### 4. Export excel and gtf for novel exons ------------------------------------
  # 4.1 Export excel table -----------------------------------------------------
# export only neoexons that passed the fold change, expression filters
neoexons = final_annot %>%
  dplyr::filter(isoform_version %in% rownames(all))
dim(neoexons)
write.xlsx(neoexons, file = ".../Neoexons_20260217.xlsx", sep="\t")

  # 4.2 Export gtf with full length isoform ------------------------------------
#.Load gtf 
gtf <- readGFF(".../Novel_isoforms_final.gtf")
dim(gtf)
class(gtf)

gtf.novel = gtf %>%
  dplyr::filter(transcript_id %in% rownames(all)) 
dim(gtf.novel)

unique(gtf.novel$gene_id)
unique(gtf.novel$transcript_id)

gtf.novel = makeGRangesFromDataFrame(gtf.novel, keep.extra.columns=TRUE,
                                     seqinfo=NULL)

rtracklayer::export(gtf.novel, ".../Neoexons_Full_length_isoforms_20260217.gtf")

  # 4.3 Export gtf with only novel exons  --------------------------------------
#.Load gtf 
Long <- readGFF(".../2024-06-27/query_GTF_1/EwNG_Long_isoform.gtf")
dim(Long)
class(Long)

Long = Long %>% dplyr::select(!cmp_ref_gene)

Truncated <- readGFF(".../2024-06-27/query_GTF_1/EwNG_Truncated_isoform.gtf")
dim(Truncated)
class(Truncated)

Neoisoform <- readGFF(".../2024-06-27/query_GTF_1/EwNG_Neoisoform.gtf")
dim(Neoisoform)
class(Neoisoform)
Neoisoform = Neoisoform %>% dplyr::select(!cmp_ref_gene)

gtf = rbind(Long, Truncated)
gtf = rbind(gtf, Neoisoform)

novel_transcript_id =  gsub("_.*", "", gtf$transcript_id)
gtf$transcript_id = novel_transcript_id
gtf.novel = gtf %>%
  dplyr::filter(transcript_id %in% rownames(all)) 
dim(gtf.novel)

unique(gtf.novel$gene_id)
unique(gtf.novel$transcript_id)

gtf.novel = makeGRangesFromDataFrame(gtf.novel, keep.extra.columns=TRUE,
                                     seqinfo=NULL)

rtracklayer::export(gtf.novel, ".../Neoexons_only_20260217.gtf")

### 5. Matrix transformation and filtering -------------------------------------
names <- rownames(all)
all = sapply(all, as.numeric)
all = as.matrix(all)
rownames(all) = names

### 6. Plot heatmap per isoform ------------------------------------------------

loadfonts()
cols <-c("#6947A4", "#FBF6F4", "red")
cols <- paletteer::paletteer_c("grDevices::Blue-Red", n=100)
mat = t(apply(all, 1, scale))
colnames(mat) <- colnames(all)

mat = t(mat)

# Extract and sort the row names
sorted_row_names <- sort(rownames(mat))

# Reorder the matrix based on sorted row names
ordered_mat <- mat[sorted_row_names, ]

colnames(ordered_mat)

# Add foldchange annotation 
fold_change_Illumina <- size
fold_change_Illumina <- fold_change_Illumina[rownames(fold_change_Illumina) %in% rownames(all), ]
fold_change_Illumina <- fold_change_Illumina$FC_Ewing_Illumina
names(fold_change_Illumina) <- rownames(all)

fold_change_Nanopore <- size
fold_change_Nanopore <- fold_change_Nanopore[rownames(fold_change_Nanopore) %in% rownames(all), ]
fold_change_Nanopore <- fold_change_Nanopore$FC_Ewing_Nanopore
names(fold_change_Nanopore) <- rownames(all)

fold_change_annotation <- rowAnnotation(
  fold_change_Illumina = log10(0.0001 + fold_change_Illumina), 
  fold_change_Nanopore = log10(0.0001 + fold_change_Nanopore),
  col = list(
    fold_change_Illumina = colorRamp2(
      c(min(log10(0.0001 + fold_change_Illumina)), 0, max(log10(0.0001 + fold_change_Illumina))), 
      c("#F1F1F1FF", "#8DA3CAFF", "#6B0077FF")),
    fold_change_Nanopore = colorRamp2(
      c(min(log10(0.0001 + fold_change_Nanopore)), 0, max(log10(0.0001 + fold_change_Nanopore))), 
      c("#F2F0F6FF", "#D68ABEFF", "#7D0112FF"))
  ),
  annotation_name_gp = gpar(fontsize = 5,  col = "black", fontfamily = "Helvetica", fontface ="plain"),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 5, fontfamily = "Helvetica", fontface = "plain"),
    labels_gp = gpar(fontsize = 5, fontfamily = "Helvetica")
  ),
  simple_anno_size = unit(0.2, "cm"), 
  width = unit(0.1, "cm")
)

cols <- c("#F4FAFEFF", "#7FABD3FF", "#273871FF")

flipped_mat <- t(ordered_mat)

h <- Heatmap(
  flipped_mat, 
  col = cols, 
  name = "Expression Z-score (TPM)",
  border_gp = gpar(col = "black", lty = 3),
  cluster_rows = FALSE,       
  cluster_columns = FALSE,
  show_row_names = FALSE,     
  show_column_names = TRUE,   
  column_names_gp = gpar(fontsize = 5, fontfamily = "Helvetica"),
  row_names_gp = gpar(fontsize = 5, fontfamily = "Helvetica"),
  heatmap_legend_param = list(          # **ONLY HERE**
    title_gp = gpar(fontsize = 5, fontfamily = "Helvetica", fontface = "plain"),
    labels_gp = gpar(fontsize = 5, fontfamily = "Helvetica"),
    legend_height = unit(2, "cm"),
    legend_width = unit(1, "cm")
  ),
  right_annotation = fold_change_annotation
)

draw(h,
     heatmap_legend_side = "right",
     annotation_legend_side = "right"
)

png(paste0(resdir, "heatmap_neogenes_per_gene_FC15_Rowmeans_1_min_exp_20.png"), 
    width = 18, height = 10, units = "cm", res = 300)
draw(h)
dev.off()

### 7. Plot heatmap per gene ---------------------------------------------------
# Calculate fold change per gene
size.2 = size
original_rownames <- rownames(size.2)
new_rownames <- sub("\\..*", "", original_rownames)
size.2$gene<- new_rownames

size.2  <- size.2 %>%
  group_by(gene) %>%
  summarize(across(everything(), mean, na.rm = TRUE)) %>%
  arrange(-FC_Ewing_Illumina)

ord <- size.2$gene

all.gene <- as.data.frame(all)

original_rownames <- rownames(all.gene)
new_rownames <- sub("\\..*", "", original_rownames)
all.gene$gene<- new_rownames

# For each gene calculate the mean expression
gene_means <- all.gene %>%
  group_by(gene) %>%
  summarize(across(everything(), mean, na.rm = TRUE))

rownames(gene_means) = gene_means$gene
names <- rownames(gene_means)

ord = ord[ord %in% rownames(gene_means)]
gene_means <- gene_means[match(ord, rownames(gene_means)), ]
rownames(gene_means) = gene_means$gene
names = rownames(gene_means)
gene_means = gene_means %>% dplyr::select(-gene)
rownames(gene_means) = names

# Retransform into matrix 
gene_means = sapply(gene_means, as.numeric)
gene_means = as.matrix(gene_means)
rownames(gene_means) = names

# Number of neogenes kept
dim(gene_means)

# Color set up 
cols <-c("#6947A4", "#FBF6F4", "red")
mat = t(apply(gene_means, 1, scale))
colnames(mat) <- colnames(gene_means)

mat = t(mat)

# Extract and sort the row names
sorted_row_names <- sort(rownames(mat))

# Reorder the matrix based on sorted row names
ordered_mat <- mat[sorted_row_names, ]

# Add foldchange annotation 
fold_change_Illumina <- size.2
fold_change_Illumina = fold_change_Illumina[fold_change_Illumina$gene %in% ord,]
fold_change_Illumina = fold_change_Illumina$FC_Ewing_Illumina
names(fold_change_Illumina) <- rownames(gene_means)

fold_change_Nanopore <- size.2
fold_change_Nanopore = fold_change_Nanopore[fold_change_Nanopore$gene %in% ord,]
fold_change_Nanopore = fold_change_Nanopore$FC_Ewing_Illumina
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
  annotation_height = unit(c(0.01, 0.01), "cm"),  # Adjust the height of the annotation bars
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
  name = "Expression Z-score (TPM)",
  border_gp = gpar(col = "black", lty = 3),
  cluster_columns = FALSE,
  show_column_names = TRUE,
  cluster_rows = FALSE,
  show_row_names = TRUE,
  show_heatmap_legend = TRUE,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  top_annotation = fold_change_annotation
)

# Draw the heatmap with the legend on the right side
draw(h, heatmap_legend_side = "right")

png(paste0(resdir, "heatmap_novel_exons_per_gene.png"), width = 30, height = 15, units = "cm", res = 300)
draw(h)
dev.off()


