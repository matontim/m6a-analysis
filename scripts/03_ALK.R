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
