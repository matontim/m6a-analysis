library(m6AnetAnalyzer)

# Loading Dataset
data_naive <- read.csv("data/naive/pooled_naive_2025_data.site_proba.csv")
head(data_naive)

# Generate summary HTML report
summarize_m6anet_output(
    data_naive,
    output_file = "naive_summary_report.html",
    output_dir = "results")

# Creating filtered dataset for downstream analysis
filtered_data_naive <- filter_m6a_sites(data_naive,
prob_modified =0.9)

# Save output
write.csv(filtered_data_naive, file.path("data/naive", "filtered_data_naive.csv"))
