[![Jacquemont's Lab Header](labheader.png)](https://www.jacquemont-lab.org/)

[Git Repository MCNV2 – Mendelian CNV Validation](https://github.com/JacquemontLab/MCNV2-Mendelian-CNV-Validation)

# MCNV2 – Mendelian CNV Validation

## Overview

MCNV2 is an R package for the validation and annotation of Mendelian CNVs. It provides tools for preprocessing CNV data, gene annotation, and visualization via a Shiny app.

---

## Installation

MCNV2 can be installed directly from GitHub:
```r
# Install devtools if not already installed
install.packages("devtools")

# Install the MCNV2 package 
devtools::install_github("JacquemontLab/MCNV2-Mendelian-CNV-Validation")
```

Or if you have a local source tarball:
```r
install.packages("/path/to/MCNV2_0.1.0.tar.gz", repos = NULL, type = "source")

---

## Requirements

### R Packages

The required packages are listed in the `DESCRIPTION` file.

---

### Python 3

MCNV2 relies on Python 3 for CNV annotation and inheritance calculation.

* Ensure that `python3` is installed and available in R's PATH.
* If you want to use a specific Python version, specify it in R:

```r
library(reticulate)
use_python(Sys.which("python3"), required = TRUE)
```

* MCNV2 also requires a few Python packages, which can be installed in a dedicated virtual environment:

```r
library(MCNV2)
MCNV2::setup_python_env(envname = "r-MCNV2")
```

This will create a virtual environment named `r-MCNV2` and install all necessary Python dependencies listed in `inst/python/requirements.txt` (e.g., `polars`).

---

### Bedtools

`bedtools` is required for certain CNV analyses.

* Make sure it is installed and in your PATH.
* To use a specific version, specify the full path in R:

```r
library(MCNV2)
MCNV2::launch(bedtools_path = Sys.which("bedtools"), results_dir = "~/projects/mcnv2_results")
```

* Alternatively, you can download a precompiled release from [MCNV2 releases](https://github.com/xxx/MCNV2-Mendelian-CNV-Validation/releases) and install from source in R:

```r
install.packages("/path/to/MCNV2_vX.XX.XXX", type = "source", repos = NULL)
```

> Note: `bedtools` is not available natively for Windows. Use a virtual machine or [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/install-win10).

---

### Operating System Support

MCNV2 has been tested on:

* macOS 15.6.1 "Sequoia"
* Linux distributions (tested on Ubuntu 22.04)

Windows users must rely on WSL or a VM for `bedtools` functionality.

---

## Example Usage

After installing dependencies and setting up Python:

```r
library(reticulate)
use_virtualenv("r-MCNV2", required = TRUE)

library(MCNV2)
launch(bedtools_path = Sys.which("bedtools"), results_dir = "~/projects/mcnv2_results")
```

This will launch the interactive Shiny app for CNV validation and annotation.


