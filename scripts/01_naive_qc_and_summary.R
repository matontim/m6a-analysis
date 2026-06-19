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
prob_modified =0.9)

# Save output
write.csv(
    filtered_data_naive, 
    file.path("data/naive", 
    "filtered_data_naive.csv")
    )

# Path to GTF annotation file 
gtf_path <- "/data/GTF_annotation"

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
write.csv(naive_m6A_with_tx_annotations, file.path("results/tables","naive_m6A_with_tx_annotations.csv"))
ggsave(file.path("results/tables","naive_tx_density.pdf"), p)


