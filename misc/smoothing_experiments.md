# Introduction

This file contains the instructions to reproduce the results of the smoothing experiments.

So far, the following methods have been tested:

- Loess (Location: `src/tox_loess.F90`)
- Nadaraya–Watson (Location: `src/knn_smoothing_nadaraya_watson.F90`)
- ANWIL (Location: `src/anwil.F90`)
- ManLe (Location: `src/manle_module.F90`)
- AManLe (Location: `src/manle_module.F90`)

At the moment, both ManLe and AManLe are implemented in the same module.

*Important: Since experiments are still ongoing, there may be allocations inside the subroutines. This will be fixed in the future.* 

Please refer to the method details in [Fortran_Coding_Guides.pdf](https://gitlab.rlp.net/a.hallab/tensor-omics/-/blob/Smoothing-Algorithm-Descriptions/misc/Tensor_Omics_Methods.pdf?ref_type=heads).

Now we're testing another approach called LoManle. This tests are independent from the previous ones so please read carefully how to run these tests.

# Project compilation

To reproduce the tests or modify the code and try it out, make sure you always compile the project.  
Apply your changes inside the `src/` folder and compile as usual using:

```bash
./build.sh
````

To run the test suite for all modules in the project, execute:

```bash
./test_runner.sh
```

# Data generation

We generate datasets to experiment with and compare the methods using a Python script. This script generates CSV files in `results/data` for different datasets with low, medium, and high noise levels. You can find the script at:

```
python/generate_smoothing_datasets.py
```

You can also find a script to plot those datasets in:

```
python/plot_datasets.py
```

# Smoothing-specific tests

To generate the smoothing-specific tests, we use a Fortran program that reads the previously generated datasets, applies smoothing using all methods, and writes CSV files for each dataset containing all relevant output from each method, with a header of the form:

```
x_original,y_original,x_loess,y_loess,x_anwil,y_anwil,x_anwil_iterative,y_anwil_iterative,x_nw,y_nw,y_nw_knn,x_manle,y_manle,x_manle_svd,y_manle_svd,x_amanle,y_amanle
```

The source file is located at:

```
test_aux/smooth_all_methods.f90
```

## Automatic compiling

### Loess, Anwil, Nadaraya Watson, Manle, Amanle

To generate smoothing-specific tests, we use a Fortran program that applies various smoothing methods (Loess, ANWIL, NW, ManLe, AManLe) to the generated datasets.

Use the script `./run_smoothing_tests.sh` to compile the project, run the smoothing across a grid of hyperparameters, and generate comparison plots. Now we are generating plots for ANWIL, Manle, AMANLE and local sigma smoothing only. Please modify the script to generate all plots or the local sigma plots separatly (uncomment the desired script).

Please consider that for the latest experiments we used complete nadaraya watson for local sigma smoothing. If you want to use knn, please modify the call `smooth_vectors_gaussian_adaptive_nw` in `anwil` module to use `smooth_type` 2.

**Usage:**
```bash
./run_smoothing_tests.sh <all|dataset_name> <k_list> <span_list> <iters_list> <kernel_list> <sigma_k_list> <method_id> <score_list> <w_rough> <w_rmse> <w_cov>

```

**Parameters:**
* `all|dataset_name`: If `all` is used, all datasets in `results/data/`  directory are going to generate the smoothed results. You can also only generate the results for one single dataset passing the name of the file.
* `k_list`: Comma-separated list of neighbors (e.g., `10,20`).
* `span_list`: Comma-separated list of spans for Loess (e.g., `0.30,0.50`).
* `iters_list`: Comma-separated list of max iterations.
* `kernel_list`: `1` for Gaussian, `2` for Tricubic.
* `sigma_k_list`: Neighbors used for local sigma smoothing in ANWIL.
* `method_id`: ID of the method to execute (0=All, 1=Loess, ..., 7=AManLe). See below.
* `score_list`: Calculate smoothing score using 1 = arithmetic mean, 2 = geometric mean
* `w_rough`: Comma-separated weights for Roughness
* `w_rmse`: Comma-separated weights for RMSE
* `w_cov`: Comma-separated weights for Coverage.


Where `method_id`:
- 0 to execute all methods
- 1 to execute loess
- 2 to execute isotropic anwil 
- 5 to execute nadaraya watson
- 6 to execute manle
- 7 to execute amanle 



The smoothed results are located under `results/data/`.

The plots are located in `results/plots/` with the corresponding parameters that you used in the call.

### Automated Manifold Learning Pipeline (`LoManLe`)

This bash script automates the compilation, execution, and visualization of the **LoManLe** (Local Manifold Learning) algorithm across multiple datasets. It handles the batch processing of CSV files and dynamically generates reports using R.

#### 1. Features

* **Automatic Compilation**: Links the Fortran source with LAPACK, BLAS, and LOESS libraries.
* **Batch Processing**: Use the `all` keyword to process every CSV in the data directory.
* **Dynamic Parameterization**: Allows testing of multiple $k$ (nearest neighbors) values in a single run.
* **Integrated Visualization**: Automatically calls R scripts to generate PDF reports after each execution.

#### 2. Usage

```bash
./run_lomanle_tests.sh <input_file|all> <k_list> [dimension] [gap] [overlap_max] [overlap_min]
```

##### Arguments

| Argument | Type | Description | Default |
| --- | --- | --- | --- |
| `input_file` | String | Path to a `.csv` file or the keyword `all`. | Mandatory |
| `k_list` | List | Comma-separated values (e.g., `10,20,30`). | Mandatory |
| `dimension` | Integer | Manifold dimension ($d=1$ for curves, $d=2$ for surfaces). | `1` |
| `gap` | Float | $G$ threshold for identifying local discontinuities. | `1.0` |
| `overlap_max` | Float | Maximum allowable overlap between local charts. | `0.30` |
| `overlap_min` | Float | Minimum connectivity required between charts. | `0.05` |

#### 3. Workflow Logic

1. **Compilation**: The script ensures `test_lomanle.f90` is compiled with the necessary headers (`-I build`) and linked against the pre-compiled objects in `build/*.o`.
2. **Execution Loop**: For every $k$ in the provided list, the script runs the Fortran binary.
3. **File Management**: It renames the generic `lomanle_output.csv` to a unique name incorporating the parameters (e.g., `data_k10_g1.0_omax0.30_lomanle.csv`).
4. **Plotting**: Executes `plot_lomanle_spheres.R`, which reads the output and generates a visual representation of the manifold patches.

#### 4. Example Command

To process a specific dataset with two different neighbor counts:

```bash
./run_lomanle_tests.sh results/data/circular_arc_noise_high.csv 30 1 3.0 0.3 0.10
```

## Manual compiling

You can compile it as follows:

```bash
./build.sh

gfortran -O2 -I build -o build/smooth_all test_aux/smooth_all_methods.f90 build/*.o -Lexternal/lib -lloess_netlib -llapack -lblas
```

To run it on a dataset, use:

```bash
./build/smooth_all results/data/circular_arc_noise_high.csv k_neighbors n_iters_max method k_neighbors_sigma kernel_type span smoothing_score_type w_r w_e w_c
```

To generate all datasets in a single run:

```bash
for f in results/data/*.csv; do
  if [[ "$f" != *"_smoothed"* ]]; then
    echo "Processing: $f"
    ./build/smooth_all "$f" k_neighbors n_iters_max method k_neighbors_sigma kernel_type span smoothing_score_type w_r w_e w_c
  fi
done
```

This skips datasets that already have the _smoothed suffix and only processes the original input datasets.

