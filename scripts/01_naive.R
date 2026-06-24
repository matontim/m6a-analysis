library(m6AnetAnalyzer)

# Loading Dataset
data_naive <- read.csv("data/naive/pooled_naive_2025_data.site_proba.csv")
head(data_naive)

# Generate summary HTML report
summarize_m6anet_output(
                        data_naive,
                        output_file = "naive_summary_report.html",
                        output_dir = "results/tables")

# Creating filtered dataset for downstream analysis
filtered_data_naive <- filter_m6a_sites(data_naive,
                                        prob_modified = 0.9)

# Save output
write.csv(filtered_data_naive, file.path("data/naive",
                                         "filtered_data_naive.csv"))

# Path to GTF annotation file
gtf_path <- "data/GTF_annotation/gencode.v45.annotation.gtf.gz"

# Step 1: Extract transcript region lengths (5'UTR, CDS, 3'UTR)
tx_annotations <- get_transcript_region_lengths(gtf_path)

# Load in annotations
data("tx_annotations")
head(tx_annotations)

# Step 2: Annotate each m6A site with its region and relative position
naive_m6A_with_tx_annotations <-
  map_relative_tx_regions_to_m6A(
    filtered_data_naive,
    tx_annotations
  )

head(naive_m6A_with_tx_annotations)

# Step 3: Visualize distribution of m6A across transcript regions
p <- plot_relative_positions(naive_m6A_with_tx_annotations, label = "naive")

# Show plot
p

# Save plot as a PNG
ggsave(file.path("results/figures", "plot_relative_positions.png"), p)

# Save output
write.csv(naive_m6A_with_tx_annotations,
          file.path("results/tables", "naive_m6A_with_tx_annotations.csv"))
ggsave(file.path("results/tables", "naive_tx_density.pdf"), p)

# Chromosomal annotation
chr <- calculate_chromosome_location(filtered_data_naive,
                                     output_csv = file.path("results/tables",
"chromosomal_distribution.csv"),
                                     output_plot = file.path("results/figures",
"naive_chr_distribution.png"))

# Show chromosome annotation
head(chr)

# Weighted modification ratio
naive_wmr <- calculate_weighted_mod_ratio(naive_m6A_with_tx_annotations,
                                          mod_ratio_column = "mod_ratio")

head(naive_wmr)

# Convert to genomic coordinates and create BED files
bed_output <- create_bed_file(
  filtered_data_naive,
  gtf_path,
  output_bed = "results/tables/m6A_naive_gli36.bed",
  output_bedgraph = "results/tables/m6A_naive_gli36.bedgraph"
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
annotated_naive_with_rbp <- annotate_m6a(
  bed_output$m6a_df_with_genoimics_coordinates,
  example_rbp_annotation,
  feature_col = "RBP_name"
)

# Show output #1: dataframe with the overlapping features and m6A row index
head(annotated_naive_with_rbp$overlap_df)