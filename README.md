# Dark_proteome_ewing_sarcoma
Available code for manuscript Oncofusion-driven transcription unlocks a multilayered neoprotein landscape in Ewing sarcoma.

# 0.Filtering_neotranscripts.R
After Nanopore and PacBio isoforms merge, this script allow for EwS-specific isoforms (neoisoforms + neogenes) selection.
This was performed per cell line.

# 1.Combine_all_models.R
The 5 cell lines processed in Nanopore and PacBio are now merged into a single file.

# 2.Recode_gtf.R
Create a gtf with neoisoforms and neogenes identified in Long-read RNAseq.

# 3.Remove_FSM_filter_neoexons.R
Remove full-splice match isoforms in comparision with Gencode and identify EwS-specific neoexons.

# 4.Create_gtf_neoisoforms_without_full_splice_match.R
Create the matching gtf with previous #3 script.

# 5.Filter_neogenes_expression_dataset.R
Keep only Neogenes that are not expressed in TCGA and GTEX.

# 6.Update_isoforms_annot.R
Re-annotate column annotations for isoforms.

# 7.GViz_TRPM4.R
Genomic visualisation using GViz package for TRPM4 gene.

# 8.Neoisoforms_analysis.R
Create upsetplots for Nanopore, PacBio and barplots with number of isoforms in each category.

# 9.Riboseq_expression.R
Neoisoforms Riboseq expression in dox-inducible cellular models.

# 10.Make_supplementary_Table_S6.R
Clean annotations for supp table S6.

# 11.Make_Neogene_neoisoforms_hg38_gtf.R
Make final neoisoforms gtf

# 12.Deeploc_neoisoforms.R
Analysis of Deeploc localization predictions for neoisoforms translons.

# 13.Interproscan.R
Analysis of Interproscan predictions for neoisoforms translons.

# 14.Do_some_statistics.R
Count isoforms.

# 15.Riboseq_expression_KASH5.R
Riboseq expression of KASH5 and other cancer testis antigens.

# 16.Riboseq_expression_TRPM4.R
Riboseq expression of TRPM4.

# 17.shRNA_Ew_NG48.R
qPCR plots of Ew_NG48 upon shRNA mediated knock-down.

# 18.CRISPRi_Ew_NG48.R
qPCR plots of Ew_NG48 upon CRISPRi targeting regions flanking FLI1-bound GGAA microsatellites.







