library(dplyr)
library(readr)

# Daten laden
results <- c("pc67_1", "pc196_1", "pc197_1", "pc202_1", "pc247_1", "pc302_1",
	"pc325_1", "pc346_1", "pc418_1", "pc483_1", "pc538_1", "pc992_1", "pc1020_1",
	"pc1088_1", "pc1375_1", "pc1446_1", "pc1643_1", "pc1675_1", "pc1885_1", "pc2317_1",
	"pc3076_1", "pc3257_1", "pc3420_1", "pc3921_1", "pc5043_1", "pc5638_1", "pc6096_1",
	"pc6500_1", "pc7658_1", "pc8012_1")

go <- readRDS("../material/protein2goa.rds")
ipr <- readRDS("../material/protein2ipr.rds")
hrd <- readRDS("../material/protein2hrd.rds")

extract_annotations <- function(annotation_data, result_list, annotation_col) {
  result_ids <- unique(result_list$ids)
  
  annotation_data %>%
    select(protein_id, {{annotation_col}}) %>%
    filter(protein_id %in% result_ids) %>%
    group_by(protein_id) %>%
    summarize(
      annotations = list(unique({{annotation_col}})),
      annotations_str = paste(unique({{annotation_col}}), collapse = "|"),
      n_annotations = n_distinct({{annotation_col}}),
      .groups = "drop"
    )
}

for (name in results) {
  message("Verarbeite: ", name)
  result <- readRDS(paste0("../results/", name, "/", name, "_output.rds"))
  
  # KORRIGIERT: Direkter Zugriff auf die IDs
  go_annotations <- extract_annotations(go, result$per_corpus$GO, GOA_label)
  ipr_annotations <- extract_annotations(ipr, result$per_corpus$IPR, description)  
  hrd_annotations <- extract_annotations(hrd, result$per_corpus$HRD, hrd)
  
  # Als CSV speichern
  write_csv(go_annotations, paste0("../results/", name, "/", name, "_go_annotations.csv"))
  write_csv(ipr_annotations, paste0("../results/", name, "/", name, "_ipr_annotations.csv")) 
  write_csv(hrd_annotations, paste0("../results/", name, "/", name, "_hrd_annotations.csv"))
  
  message("  GO: ", nrow(go_annotations), " Zeilen")
  message("  IPR: ", nrow(ipr_annotations), " Zeilen")
  message("  HRD: ", nrow(hrd_annotations), " Zeilen")
}
