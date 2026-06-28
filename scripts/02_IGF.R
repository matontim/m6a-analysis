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