source("r/load_tensor_omics.R")
source("r/test_helpers.R")

# N=4, 2 ensembles: {1,2,3} (seed=1) and {2,3,4} (seed=4), overlapping on {2,3} -- Overlap
# Coefficient 2/3 -- merged into one super-ensemble. Point 1 is also flagged low-confidence.
# Matches mod_test_shape_truthful_clustering_json.R's own fixture.
fixture <- function() {
  seed_selection_mask <- c(TRUE, FALSE, FALSE, TRUE)
  ensemble_masks <- matrix(FALSE, nrow = 4, ncol = 2)
  ensemble_masks[1:3, 1] <- TRUE
  ensemble_masks[2:4, 2] <- TRUE
  ensemble_low_confidence_masks <- matrix(FALSE, nrow = 4, ncol = 2)
  ensemble_low_confidence_masks[1, 1] <- TRUE
  super_ensembles <- matrix(c(1, 2, 0, 0), nrow = 2, ncol = 2)
  list(seed_selection_mask = seed_selection_mask, ensemble_masks = ensemble_masks,
       ensemble_low_confidence_masks = ensemble_low_confidence_masks, super_ensembles = super_ensembles)
}

test_points_csv_membership <- function() {
  fx <- fixture()
  filename <- "test_stc_points_r.csv"
  serialize_stc_points_as_csv(filename, 1L, fx$seed_selection_mask, fx$ensemble_masks,
                               fx$ensemble_low_confidence_masks, fx$super_ensembles)

  lines <- readLines(filename, warn = FALSE)
  file.remove(filename)

  assert_true(lines[1] == "row,ensembles,super_ensembles,low_confidence_ensembles,seed_of")
  assert_true(lines[2] == '1,"1","1","1","1"')
  assert_true(lines[3] == '2,"1,2","1","",""')
  assert_true(lines[4] == '3,"1,2","1","",""')
  assert_true(lines[5] == '4,"2","1","","2"')
}

test_points_csv_zero_ensembles <- function() {
  seed_selection_mask <- c(FALSE, FALSE)
  ensemble_masks <- matrix(FALSE, nrow = 2, ncol = 0)
  ensemble_low_confidence_masks <- matrix(FALSE, nrow = 2, ncol = 0)
  super_ensembles <- matrix(0L, nrow = 2, ncol = 0)

  filename <- "test_stc_points_zero_r.csv"
  serialize_stc_points_as_csv(filename, 0L, seed_selection_mask, ensemble_masks,
                               ensemble_low_confidence_masks, super_ensembles)

  lines <- readLines(filename, warn = FALSE)
  file.remove(filename)

  assert_true(lines[2] == '1,"","","",""')
}

test_overlap_csv_pairwise_coefficient <- function() {
  fx <- fixture()
  filename <- "test_stc_overlap_r.csv"
  serialize_stc_ensemble_overlap_as_csv(filename, fx$ensemble_masks)

  lines <- readLines(filename, warn = FALSE)
  file.remove(filename)

  assert_true(lines[1] == "ensemble_a,ensemble_b,overlap_coefficient")
  assert_true(lines[2] == "1,2,6.6666666666666663E-001")
}

test_overlap_csv_no_intersecting_pairs <- function() {
  ensemble_masks <- matrix(FALSE, nrow = 4, ncol = 2)
  ensemble_masks[1:2, 1] <- TRUE
  ensemble_masks[3:4, 2] <- TRUE

  filename <- "test_stc_overlap_disjoint_r.csv"
  serialize_stc_ensemble_overlap_as_csv(filename, ensemble_masks)

  lines <- readLines(filename, warn = FALSE)
  file.remove(filename)

  non_empty <- lines[nzchar(lines)]
  assert_true(length(non_empty) == 1)
  assert_true(non_empty[1] == "ensemble_a,ensemble_b,overlap_coefficient")
}

test_super_ensembles_tsv_gene_family_format <- function() {
  super_ensembles <- matrix(c(1, 2, 0), nrow = 3, ncol = 1)

  filename <- "test_stc_super_ensembles_r.tsv"
  serialize_stc_super_ensembles_as_tsv(filename, super_ensembles)

  lines <- readLines(filename, warn = FALSE)
  file.remove(filename)

  assert_true(lines[1] == "1\t1,2")
}

test_super_ensembles_tsv_empty_when_zero_groups <- function() {
  super_ensembles <- matrix(0L, nrow = 2, ncol = 0)

  filename <- "test_stc_super_ensembles_empty_r.tsv"
  serialize_stc_super_ensembles_as_tsv(filename, super_ensembles)

  lines <- readLines(filename, warn = FALSE)
  file.remove(filename)

  non_empty <- lines[nzchar(lines)]
  assert_true(length(non_empty) == 0)
}

run_all_tests()
