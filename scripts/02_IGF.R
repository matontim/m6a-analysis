library(m6AnetAnalyzer)

# Loading Dataset
data_igf <- read.csv("data/igf2bp2_kd/pooled_igf_2025_m6a.csv")
head(data_igf)

# Generate summary HTML report
summarize_m6anet_output(
                        data_igf,
                        output_file = "igf_summary_report.html",
                        output_dir = "results/tables")

# Creating filtered dataset for downstream analysis
filtered_data_igf <- filter_m6a_sites(data_igf,
                                      prob_modified = 0.9)

# Save output
write.csv(filtered_data_igf, file.path("data/igf2bp2_kd",
                                       "filtered_data_igf.csv"))

# Path to GTF annotation file
gtf_path <- "data/GTF_annotation/gencode.v45.annotation.gtf.gz"

# Step 1: Extract transcript region lengths (5'UTR, CDS, 3'UTR)
tx_annotations <- get_transcript_region_lengths(gtf_path)

# Load in annotations
data("tx_annotations")
head(tx_annotations)

# Step 2: Annotate each m6A site with its region and relative position
igf_m6A_with_tx_annotations <-
  map_relative_tx_regions_to_m6A(
    filtered_data_igf,
    tx_annotations
  )
head(igf_m6A_with_tx_annotations)

# Step 3: Visualize distribution of m6A across transcript regions
p <- plot_relative_positions(igf_m6A_with_tx_annotations, label = "igf2bp2_kd")

# Show plot
p

# Save plot as a PNG
ggsave(file.path("results/figures", "plot_relative_positions_igf.png"), p)

# Save output
write.csv(igf_m6A_with_tx_annotations,
          file.path("results/tables", "igf_m6A_with_tx_annotations.csv"))
ggsave(file.path("results/tables", "igf_tx_density.pdf"), p)

# Chromosomal annotation
chr <- calculate_chromosome_location(filtered_data_igf,
                                     output_csv = file.path("results/tables",
"chromosomal_distribution_igf.csv"),
                                     output_plot = file.path("results/figures",
"igf_chr_distribution.png"))

# Show chromosome annotation
head(chr)

# Weighted modification ratio
igf_wmr <- calculate_weighted_mod_ratio(igf_m6A_with_tx_annotations,
                                        mod_ratio_column = "mod_ratio")

head(igf_wmr)

# Convert to genomic coordinates and create BED files
bed_output_igf <- create_bed_file(
  filtered_data_igf,
  gtf_path,
  output_bed = "results/tables/m6A_igf_gli36.bed",
  output_bedgraph = "results/tables/m6A_igf_gli36.bedgraph"
)

# Load example RBP annotation
data("example_rbp_annotation")

# Define columns
colnames(example_rbp_annotation) <- c(
  "chr", "start", "end", "peak_id", "strand",
  "RBP_name", "experiment_method", "sample",
  "accession_of_raw_data", "conf_score"
)

# Show first few lines
head(example_rbp_annotation)

# Step 3: Run annotation function
# Note: annotate_m6a() expects columns chr, start, end, strand
# bed_output$bed_df (BED6 format) has these column names, whereas
# bed_output$m6a_df_with_genomic_coordinates uses 'genomic_pos' instead
# of 'start'/'end' and fails the function's column check.
annotated_igf_with_rbp <- annotate_m6a(
  bed_output_igf$bed_df,
  example_rbp_annotation,
  feature_col = "RBP_name"
)

# Show output #1: dataframe with the overlapping features and m6A row index
head(annotated_igf_with_rbp$overlap_df)

# Show output #2: m6a dataframe with the summarized overlaping features
annotated_m6A_df_igf <- annotated_igf_with_rbp$annotated_m6A_df_igf
head(annotated_m6A_df_igf[!is.na(annotated_m6A_df_igf$Overlapping_Features), ])

# Using built in example external datasets
# Note: example built in datasets were already demonstrated in 01_naive.R
# File paths
rbp_file <- system.file("extdata", "RBP_POSTAR_human_subset.bed",
package = "m6AnetAnalyzer")
snp_file <- system.file("extdata", "dbSNP_human_subset.bed",
package = "m6AnetAnalyzer")

# Load data
rbp_df <- readr::read_tsv(rbp_file, col_names = F)
snp_df <- readr::read_tsv(snp_file, col_names = F)

# set columns
colnames(rbp_df) <- c(
  "chr",                   # Chromosome of the RBP binding site
  "start",                 # Start position (0-based)
  "end",                   # End position
  "peak_id",               # Unique identifier for the peak
  "strand",                # Strand of the peak (+ or -)
  "RBP_name",              # Name of the RNA-binding protein
  "experiment_method",     # Method used to detect the peak (e.g., CLIP-seq, PIP-seq)
  "sample",                # Sample or condition in which the peak was observed
  "accession_of_raw_data", # Accession number for the raw experiment data
  "conf_score"             # Confidence score or enrichment measure for the peak
)

colnames(snp_df) <- c(
  "chr",                  # Chromosome of the SNP
  "start",                # Start position (0-based)
  "end",                  # End position (1-based; SNP occupies 1 bp)
  "snp_id",               # Unique SNP identifier
  "variation",            # Allele variation (e.g., A>G)
  "strand",               # Strand information (usually "*" for SNPs)
  "clinical_significance",# Clinical annotation of the SNP (e.g., benign, pathogenic)
  "function_class",       # Functional classification (e.g., exonic, intronic)
  "gene"                  # Associated gene symbol or ID
)

head(rbp_df)

head(snp_df)