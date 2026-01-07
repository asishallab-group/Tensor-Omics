# Introduction

This file contains the instructions to reproduce the results of the smoothing experiments.

So far, the following methods have been tested:

- Loess (Location: `src/tox_loess.F90`)
- Nadaraya–Watson (Location: `src/knn_smoothing_nadaraya_watson.F90`)
- ANWIL (Location: `src/anwil.F90`)
- ManLe (Location: `src/manle_module.F90`)
- AManLe (Location: `src/manle_module.F90`)

At the moment, both ManLe and AManLe are implemented in the same module.

Please refer to the method details in [Fortran_Coding_Guides.pdf](https://gitlab.rlp.net/a.hallab/tensor-omics/-/blob/Smoothing-Algorithm-Descriptions/misc/Tensor_Omics_Methods.pdf?ref_type=heads).

# Project compilation

To reproduce the tests or modify the code and try it out, make sure you always compile the project.  
Apply your changes inside the `src/` folder and compile as usual using:

```bash
./build.sh
````

To run the test suite for all modules in the project, execute:

```bash
./run_test.sh
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
x_original,y_original,x_loess,y_loess,x_anwil,y_anwil,x_anwil_mode1,y_anwil_mode1,x_anwil_mode2,y_anwil_mode2,x_anwil_iterative,y_anwil_iterative,x_nw,y_nw,x_manle,y_manle,x_manle_svd,y_manle_svd,x_amanle,y_amanle
```

The source file is located at:

```
test_aux/smooth_all_methods.f90
```

## Automatic compiling

Use the script `./run_smoothing_tests.sh` to compile the project, compile the program, run the tests and plot in a single step.

Usage: 

```
./run_smoothing_tests.sh <all|dataset_name> <k_neighbors> <n_iters_max>
```

If `all` is used, all datasets are going to generate the smoothed results. You can also only generate the results for one single dataset passing the name of the file.

The smoothed results are located under `results/data/`.

The plots are located in `results/plots/all_smoothed_subplots.pdf`.

## Manual compiling

You can compile it as follows:

```bash
gfortran -O2 -I build -o build/smooth_all test_aux/smooth_all_methods.f90 build/*.o -L/usr/local/lib -llapack -lblas
```

To run it on a dataset, use:

```bash
./build/smooth_all results/data/circular_arc_noise_high.csv <k_neighbors> <max_iterations>
```

This will generate a CSV file called:

```
results/data/circular_arc_noise_high_smoothed_kXX_iterYY.csv
```

with the header shown above.

To generate all datasets in a single run:

```bash
for f in results/data/*.csv; do
  if [[ "$f" != *"_smoothed"* ]]; then
    echo "Processing: $f"
    ./build/smooth_all "$f"
  fi
done
```

This skips datasets that already have the _smoothed suffix and only processes the original input datasets.