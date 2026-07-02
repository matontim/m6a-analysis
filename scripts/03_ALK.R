library(m6AnetAnalyzer)

# Loading Dataset
data_alk <- read.csv("data/alkbh5_kd/pooled_alk_2025_data.site_proba.csv")
head(data_alk)

# Generate summary HTML report
summarize_m6anet_output(
                        data_alk,
                        output_file = "alk_summary_report.html",
                        output_dir = "results/tables")

# Creating filtered dataset for downstream analysis
filtered_data_alk <- filter_m6a_sites(data_alk,
                                      prob_modified = 0.9)

# Save output
write.csv(filtered_data_alk, file.path("data/alkbh5_kd",
                                       "filtered_data_alk.csv"))

# Path to GTF annotation file
gtf_path <- "data/GTF_annotation/gencode.v45.annotation.gtf.gz"

# Step 1: Extract transcript region lengths (5'UTR, CDS, 3'UTR)
tx_annotations <- get_transcript_region_lengths(gtf_path)

# Load in annotations
data("tx_annotations")
head(tx_annotations)

# Step 2: Annotate each m6A site with its region and relative position
alk_m6A_with_tx_annotations <-
  map_relative_tx_regions_to_m6A(
    filtered_data_alk,
    tx_annotations
  )
head(alk_m6A_with_tx_annotations)

# Step 3: Visualize distribution of m6A across transcript regions
p <- plot_relative_positions(alk_m6A_with_tx_annotations, label = "alkbh5_kd")

# Show plot
p

# Save plot as a PNG
ggsave(file.path("results/figures", "plot_relative_positions_alk.png"), p)

# Save output
write.csv(alk_m6A_with_tx_annotations,
          file.path("results/tables", "alk_m6A_with_tx_annotations.csv"))
ggsave(file.path("results/tables", "alk_tx_density.pdf"), p)

# Chromosomal annotation
chr <- calculate_chromosome_location(filtered_data_alk,
                                     output_csv = file.path("results/tables",
"chromosomal_distribution_alk.csv"),
                                     output_plot = file.path("results/figures",
"alk_chr_distribution.png"))

# Show chromosome annotation
head(chr)

# Weighted modification ratio
alk_wmr <- calculate_weighted_mod_ratio(alk_m6A_with_tx_annotations,
                                        mod_ratio_column = "mod_ratio")

head(alk_wmr)

# Convert to genomic coordinates and create BED files
bed_output_alk <- create_bed_file(
  filtered_data_alk,
  gtf_path,
  output_bed = "results/tables/m6A_alk_gli36.bed",
  output_bedgraph = "results/tables/m6A_alk_gli36.bedgraph"
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
annotated_alk_with_rbp <- annotate_m6a(
  bed_output_alk$bed_df,
  example_rbp_annotation,
  feature_col = "RBP_name"
)

# Show output #1: dataframe with the overlapping features and m6A row index
head(annotated_alk_with_rbp$overlap_df)

# Show output #2: m6a dataframe with the summarized overlaping features
annotated_m6A_df_alk <- annotated_alk_with_rbp$annotated_m6A_df_alk
head(annotated_m6A_df_alk[!is.na(annotated_m6A_df_alk$Overlapping_Features), ])
