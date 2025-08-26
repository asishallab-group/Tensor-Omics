# Version 1.1.3
# Author: Aaron Schroeder

library(dplyr)
library(readr)

# List of results
results <- c("pc67_1", "pc196_1", "pc197_1", "pc202_1", "pc247_1", "pc302_1",
	"pc325_1", "pc346_1", "pc418_1", "pc483_1", "pc538_1", "pc992_1", "pc1020_1",
	"pc1088_1", "pc1375_1", "pc1446_1", "pc1643_1", "pc1675_1", "pc1885_1", "pc2317_1",
	"pc3076_1", "pc3257_1", "pc3420_1", "pc3921_1", "pc5043_1", "pc5638_1", "pc6096_1",
	"pc6500_1", "pc7658_1", "pc8012_1")

go <- readRDS("../material/protein2goa.rds")
ipr <- readRDS("../material/protein2ipr.rds")
hrd <- readRDS("../material/protein2hrd.rds")

extract_annotations <- function(annotation_data, result_list, annotation_col) {
  # keep order from k-d tree
  result_ids <- result_list$ids
  result_similarities <- result_list$similarities
  
  similarity_df <- data.frame(
    protein_id = factor(result_ids, levels = result_ids),  # ← Levels setzen!
    cosine_similarity = result_similarities
  )
  
  annotation_data %>%
    select(protein_id, {{annotation_col}}) %>%
    filter(protein_id %in% result_ids) %>%

    mutate(protein_id = factor(protein_id, levels = result_ids)) %>%
    group_by(protein_id) %>%
    summarize(
      annotations_str = paste(unique({{annotation_col}}), collapse = "|"),
      n_annotations = n_distinct({{annotation_col}}),
      .groups = "drop"
    ) %>%

    left_join(similarity_df, by = "protein_id") %>%

    select(protein_id, cosine_similarity, everything())
}

for (name in results) {
  message("Processing: ", name)
  result <- readRDS(paste0("../results_no_idf/", name, "/", name, "_output.rds"))
  

  go_annotations <- extract_annotations(go, result$per_corpus$GO, GOA_label)
  ipr_annotations <- extract_annotations(ipr, result$per_corpus$IPR, description)  
  hrd_annotations <- extract_annotations(hrd, result$per_corpus$HRD, hrd)
  
  # Save as csv
  write_csv(go_annotations, paste0("../results_no_idf/", name, "/", name, "_no_idf_go_annotations.csv"))
  write_csv(ipr_annotations, paste0("../results_no_idf/", name, "/", name, "_no_idf_ipr_annotations.csv")) 
  write_csv(hrd_annotations, paste0("../results_no_idf/", name, "/", name, "_no_idf_hrd_annotations.csv"))
  
  message("  GO: ", nrow(go_annotations), " reference proteins")
  message("  IPR: ", nrow(ipr_annotations), " reference proteins")
  message("  HRD: ", nrow(hrd_annotations), " reference proteins")
}
