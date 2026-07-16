# The script aims to create a visualizaion of TRPM4 gene 
# THe same principal was used to generate visualization for other genes TNS1, TRPM4, KASH5, NKX2-2

Sys.unsetenv("http_proxy")
Sys.unsetenv("https_proxy")
Sys.unsetenv("HTTP_PROXY")
Sys.unsetenv("HTTPS_PROXY")

### 0. Load libraries and paths ------------------------------------------------
# Load libraries
library(txdbmaker)
library(Gviz)
library(GenomicInteractions)
library(rtracklayer)
library(dplyr)

### 1. Create novel transcripts annotation -------------------------------------
# Path to the GTF file
neotranscripts <- ".../20260415_neoisoforms_hg38.gtf"

# gene coordinates in hg38
my_gene <- GRanges(
  seqnames="chr19",
  IRanges(
    start=49155884,
    end=49215008
  )
)

# Load GTF file as a TxDb object
TxDb <- txdbmaker::makeTxDbFromGFF(neotranscripts)

# Gviz gene track
neotranscriptsTrack <- GeneRegionTrack(
  TxDb,
  stacking = "squish",
  # Print gene symbols
  transcriptAnnotation="transcipt",
  # Name of the track
  name = "Genes",
  # Gene name in italic
  fontface.group="italic",
  # Remove borders around exons
  col = 0,
  # Color of the exons
  fill = "#585858",
  # Apply the exon color to the thin line in introns
  col.line = NULL,
  # Color of the gene names
  fontcolor.group= "#333333",
  # Font size
  fontsize.group=18
)

# Plot the region containing the NKX2-2 gene
plotTracks(
  neotranscriptsTrack,
  from = start(my_gene),
  to = end(my_gene),
  chromosome = as.character(seqnames(my_gene)),
  # Remove the grey background and the white borders of the track name
  background.title = "transparent",
  col.border.title="transparent",
  # Track name color
  col.title = "#333333"
)

### 2. Create Refseq annotation ------------------------------------------------
gencode <- ".../gencode.v48.annotation.gtf"

# Load GTF file as a TxDb object
TxDb <- txdbmaker::makeTxDbFromGFF(gencode)

# Load the GTF file with rtracklayer
genome_gtf <- rtracklayer::import(gencode)

# Extract gene_id and gene_symbol and remove duplicates
gene2symbol <- unique(mcols(genome_gtf)[,c("gene_id","gene_name")])

# Define gene_id as rownames
rownames(gene2symbol) <- gene2symbol$gene_id

# Gviz gene track
gencodeTrack <- GeneRegionTrack(
  TxDb,
  # collapseTranscripts = "meta",
  # Print gene symbols
  transcriptAnnotation="symbol",
  stacking = "squish",
  # Name of the track
  name = "Genes",
  # Gene name in italic
  fontface.group="italic",
  # Remove borders around exons
  col = 0,
  # Color of the exons
  fill = "#585858",
  # Apply the exon color to the thin line in introns
  col.line = NULL,
  # Color of the gene names
  fontcolor.group= "#333333",
  # Font size
  fontsize.group=18
)
ranges(gencodeTrack)$symbol <- gene2symbol[ranges(gencodeTrack)$gene, "gene_name"]

# Plot the region containing the Tram1 gene
plotTracks(
  gencodeTrack,
  from = start(my_gene),
  to = end(my_gene),
  chromosome = as.character(seqnames(my_gene)),
  # Remove the grey background and the white borders of the track name
  background.title = "transparent",
  col.border.title="transparent",
  # Track name color
  col.title = "#333333"
)

### 3. Add bigwig files --------------------------------------------------------
# path to the bigWig files
bw_folder <- ".../Bigwig/"

setwd(".../Bigwig/")
gr_orig <- import("EWSFW.bw")

# Riboseq bigwig is not annotated with ucsc, then add "chr" prefix to seqnames
seqlevels(gr_orig) <- paste0("chr", seqlevels(gr_orig))

bw_files <- c("EWSFW-ucsc.bw",
              "FLI1.bw",
              "GGAAm.bw",
              "H3K27ac.bw",
              "H3K4me3.bw",
              "Nanopore.bw", 
              "PacBio.bw",
              "ASP14D0PacBiohg38.bw",
              "ASP14D7PacBiohg38.bw")

# Extract the condition from the bw file name, i.e. the string before the first "_"
conditions <- sapply(strsplit(bw_files, "_"), `[`, 1)

colors <- c("EWSFW-ucsc.bw" = "#94d574",
            "FLI1.bw"       = "#ef3f25",
            "GGAAm.bw"      = "#fdb049",
            "H3K27ac.bw"    = "#6bb1d1", 
            "H3K4me3.bw"    = "#fec2ec",
            "Nanopore.bw"   = "#9d8af7", 
            "PacBio.bw"     = "#5f4780", 
            "ASP14D0PacBiohg38.bw" = "#80CBC4",
            "ASP14D7PacBiohg38.bw" = "#004D40")
names(colors) <- unique(conditions)

# Import the signal around the gene of interest from the bigwig files
import_bw_file <- function(folder, file, locus){
  old_wd <- getwd()
  setwd(folder)
  on.exit(setwd(old_wd))
  
  gr <- rtracklayer::import(file, which = locus)
  gr$condition <- sub("_.*", "", file)
  return(gr)
}

# For each bigWig file, make a GRanges object of the region around our gene of interest
gr_list <- lapply(bw_files, function(file){
  import_bw_file(bw_folder, file, my_gene)
})
names(gr_list) <- conditions

# Generate the bigWig signal track
# fixed = TRUE  → ylim is locked globally (GGAAm use case)
# fixed = FALSE → ylim adapts to the visible region (all other tracks)
bigwig_tracks <- function(gr, locus, colors, scale_max = NULL, fixed = FALSE){
  options(ucscChromosomeNames = FALSE)
  condition <- unique(gr$condition)
  color     <- colors[condition]
  
  if(is.null(gr$score)) stop("GRanges object must have a 'score' column")
  
  max_val <- if(!is.null(scale_max)) scale_max else max(gr$score, na.rm = TRUE)
  max_val <- ceiling(max_val / 10) * 10
  
  gTrack <- DataTrack(
    range          = gr,
    type           = "hist",
    aggregation    = "max",
    col.histogram  = 0,
    fill.histogram = color,
    baseline       = 0,
    col.baseline   = color,
    lwd.baseline   = 1,
    window         = -1,        # -1 = no binning, show raw signal
    chromosome     = as.character(seqnames(locus)),
    name           = unique(gr$condition),
    ylim           = c(0, max_val),
    yTicksAt       = seq(0, max_val, length.out = 5),
    size           = 1,
    transformation = if(fixed) NULL else identity,
    scale          = if(fixed) FALSE else TRUE
  )
  return(gTrack)
}

bw_track_list <- lapply(names(gr_list), function(name){
  if(name == "GGAAm.bw") {
    bigwig_tracks(gr_list[[name]], my_gene, colors, scale_max = 10, fixed = TRUE)
  }  else {
    bigwig_tracks(gr_list[[name]], my_gene, colors, fixed = FALSE)
  }
})
names(bw_track_list) <- names(gr_list)

# Define your desired order
desired_order <- c("FLI1.bw", "GGAAm.bw", "H3K4me3.bw", "H3K27ac.bw",
                   "Nanopore.bw", "PacBio.bw", "EWSFW-ucsc.bw",
                   "ASP14D0PacBiohg38.bw", "ASP14D7PacBiohg38.bw")

bw_track_list <- bw_track_list[desired_order]

# Plot -------------------------------------------------------------------------
plotTracks(
  c(bw_track_list, gencodeTrack, neotranscriptsTrack),
  from               = start(my_gene),
  to                 = end(my_gene),
  chromosome         = as.character(seqnames(my_gene)),
  background.title   = "transparent",
  col.border.title   = "transparent",
  col.title          = "#333333",
  col.axis           = "#333333",
  fontsize.group     = 11,
  cex.title          = 0.3,
  cex.axis           = 0.54
)


pdf(".../TRPM4_gviz.pdf", width=5, height=4)
plotTracks(
  c(bw_track_list, gencodeTrack, neotranscriptsTrack),
  from               = start(my_gene),
  to                 = end(my_gene),
  chromosome         = as.character(seqnames(my_gene)),
  background.title   = "transparent",
  col.border.title   = "transparent",
  col.title          = "#333333",
  col.axis           = "#333333",
  fontsize.group     = 11,
  cex.title          = 0.3,
  cex.axis           = 0.54
)
dev.off()

