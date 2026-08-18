source("rcpp/tensoromics_functions.R")
source("rcpp/test_helpers.R")

test_flyer_serialization <- function() {

    # -------------------------------
    # Parameters
    # -------------------------------
    n_tissues  <- 2
    n_families <- 4
    n_genes    <- 5

    filename <- "test_flyer_serialization.json"

    tissues <- c("Adipose", "Thyroid")
    family_ids <- c("EMPTYFAM1", "OG0000000", "OG0000001", "EMPTYFAM2")

    gene_ids <- c(
        "NP_001243379.1",
        "UNASSIGNED.1",
        "XP_038381480.1",
        "NP_001243386.1",
        "XP_038421312.1"
    )

    gene_species <- c(
        "Canis_lupus_protein.1",
        "Whatever_protein",
        "Canis_lupus_protein.2",
        "Canis_lupus_protein.3",
        "Canis_lupus_protein.4"
    )

    gene_types <- c("ortholog", "ortholog", "paralog", "ortholog", "paralog")

    gene_to_fam <- c(2L, 0L, 3L, 2L, 3L)
    sorted_gene_to_fam_perm <- c(2L, 1L, 4L, 3L, 5L)
    gene_outliers <- c(TRUE, FALSE, FALSE, FALSE, TRUE)

    # -------------------------------
    # Gene matrix (Fortran order → column‑major in R)
    # -------------------------------
    genes <- matrix(
        c(
            0.0541530138142222, 0.0041979991981664,
            3.1415926535897932, 2.7182818284590452,
            1.1156201021350000, 0.3080014399232190,
            0.0060563875402208, 0.1473776848574070,
            0.0060563875402208, 0.0041979991981664
        ),
        nrow = 2, ncol = 5, byrow = FALSE
    )

    # -------------------------------
    # Centroids
    # -------------------------------
    centroids <- matrix(
        c(
            0.0, 0.0,
            0.03010470067722, 0.07578784202779,
            0.56083824483761, 0.15609971956069,
            0.0, 0.0
        ),
        nrow = 2, ncol = 4, byrow = FALSE
    )

    # -------------------------------
    # Call the R wrapper
    # -------------------------------
    serialize_tox_data_as_flyer_json(
        filename,
        tissues,
        family_ids,
        gene_ids,
        gene_species,
        gene_types,
        centroids,
        genes,
        gene_to_fam,
        sorted_gene_to_fam_perm,
        gene_outliers
    )

    assert_true(file.exists(filename), "No json file created.")

    # -------------------------------
    # Load JSON and check fragment
    # -------------------------------
    text <- readLines(filename)

    expected_fragment <- paste0(
        '{',
            '"tissues":["Adipose","Thyroid"],',
            '"families":[',
                '{"family":"EMPTYFAM1","gene_indices":[],"centroid":[0.0000000000000000E+000,0.0000000000000000E+000]},',
                '{"family":"OG0000000","gene_indices":[1,4],"centroid":[3.0104700677220000E-002,7.5787842027790001E-002]},',
                '{"family":"OG0000001","gene_indices":[3,5],"centroid":[5.6083824483761002E-001,1.5609971956068999E-001]},',
                '{"family":"EMPTYFAM2","gene_indices":[],"centroid":[0.0000000000000000E+000,0.0000000000000000E+000]}',
            '],',
            '"genes":[',
                '{"coordinates":[5.4153013814222203E-002,4.1979991981664000E-003],"id":"NP_001243379.1","family":"OG0000000","species":"Canis_lupus_protein.1","is_outlier":true,"type":"ortholog"},',
                '{"coordinates":[3.1415926535897931E+000,2.7182818284590451E+000],"id":"UNASSIGNED.1","family":null,"species":"Whatever_protein","is_outlier":false,"type":"ortholog"},',
                '{"coordinates":[1.1156201021350001E+000,3.0800143992321899E-001],"id":"XP_038381480.1","family":"OG0000001","species":"Canis_lupus_protein.2","is_outlier":false,"type":"paralog"},',
                '{"coordinates":[6.0563875402208003E-003,1.4737768485740699E-001],"id":"NP_001243386.1","family":"OG0000000","species":"Canis_lupus_protein.3","is_outlier":false,"type":"ortholog"},',
                '{"coordinates":[6.0563875402208003E-003,4.1979991981664000E-003],"id":"XP_038421312.1","family":"OG0000001","species":"Canis_lupus_protein.4","is_outlier":true,"type":"paralog"}',
            ']',
        '}'
    )

    assert_true(identical(expected_fragment, text), "JSON serialization doesn't match the expected fragment.")

    file.remove(filename)
}

run_all_tests()
