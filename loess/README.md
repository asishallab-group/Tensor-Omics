# LOESS Smoothing

A Fortran implementation of the LOESS (Locally Estimated Scatterplot Smoothing) algorithm for non-parametric regression and data smoothing.

## Features

- LOESS smoothing for 1D and vector-valued data
- Outlier suppression and fallback behavior
- Well-documented Fortran source code
- Automated documentation generation with [FORD](https://github.com/Fortran-FOSS-Programmers/ford)
- Includes a test program with several validation cases

## Directory Structure

```
.
├── bin/         # Compiled files
├── src/         # Fortran source files
│   └── smoothing.f90
├── doc/         # Generated documentation (after running FORD)
├── ford.yml     # FORD configuration file
└── README.md    # Project overview
```

## Usage

### 1. Compile and run the test program

```sh
cd src
./compile_files.sh
./../bin/smoothing > output.txt
```

This will run several test cases:
- **Test Case 1:** Constant input
- **Test Case 2:** Linear trend recovery
- **Test Case 3:** Outlier suppression
- **Test Case 4:** Sparse fallback behavior
- **Test Case 5:** Vector field smoothing

Each test prints an error value and whether it passed within the expected tolerance.

### 2. Generate documentation with FORD

```sh
cd ..
ford ford.yml
```
The documentation will be generated in the `doc/` directory.
This does not work currently for an unknown reason.

## Troubleshooting

- If FORD reports "No source files with appropriate extension found", ensure that:
  - The `src/` directory exists and contains `smoothing.f90`
  - The file extension is exactly `.f90`
  - The `src_dir` in `ford.yml` is set to `src/`
- If you add more source files, place them in the `src/` directory.

## Requirements

- Fortran compiler (e.g., gfortran)
- [FORD](https://github.com/Fortran-FOSS-Programmers/ford) for documentation (optional)

## Author

Aaron Schroeder

## Version 1.0
