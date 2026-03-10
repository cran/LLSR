[![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/LLSR)](https://cran.r-project.org/package=LLSR) [![metacran downloads](https://cranlogs.r-pkg.org/badges/last-week/LLSR)](https://cran.r-project.org/package=LLSR) [![metacran downloads](https://www.r-pkg.org/badges/version/LLSR)](https://cran.r-project.org/package=LLSR) [![metacran downloads](https://www.r-pkg.org/badges/ago/LLSR)](https://cran.r-project.org/package=LLSR)

[![R-CMD-check](https://github.com/diegofcoelho/LLSR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/diegofcoelho/LLSR/actions/workflows/R-CMD-check.yaml) [![lint](https://github.com/diegofcoelho/LLSR/actions/workflows/lint.yaml/badge.svg)](https://github.com/diegofcoelho/LLSR/actions/workflows/lint.yaml) [![test-coverage](https://github.com/diegofcoelho/LLSR/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/diegofcoelho/LLSR/actions/workflows/test-coverage.yaml)

# LLSR
LLSR: Data Analysis of Liquid-Liquid Systems using R
#
Originally designed to characterise Aqueous Two Phase Systems, LLSR provide a simple way to analyse experimental data and obtain phase diagram parameters, among other properties, systematically. The package will include (every other update) new functions in order to comprise useful tools in liquid-liquid extraction research.


Supplementary Materials for "LLSR—An R Package for Data Analysis of Aqueous Two-Phase Systems", Journal of Chemical Engineering Data (2019).

Diego F. Coêlho, Jonathan G. Huddleston, Louise L. Tundisi, Pedro V. O. Menezes, Carla C. S. Porto, Priscila G. Mazzola, Roberto R. Souza and Elias B. Tambourgi

doi: [10.1021/acs.jced.9b00317](https://doi.org/10.1021/acs.jced.9b00317)

## Development

Dependencies are managed via the package `DESCRIPTION` (no renv). From the package root:

- **Run tests (installs deps if needed):** `Rscript run_tests.R` or `R -e "source('run_tests.R')"`
- **Or manually:** `remotes::install_deps()` then `testthat::test_local(".")` or `devtools::test()`

The test suite’s `tests/testthat/setup.R` also installs any missing Imports/Suggests when tests run via `R CMD check`.

