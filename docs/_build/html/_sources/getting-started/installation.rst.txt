Installation
============

MCNV2 is an **R package** that provides a **Shiny application** and **Python-backed tools** for Mendelian CNV validation and annotation.

Python is used for CNV annotation and Mendelian classification. bedtools is required for genomic overlap calculations.

.. note::

   MCNV2 has been tested on **macOS** and **Linux**.
   On Windows, we recommend using **WSL** (Windows Subsystem for Linux).

1) Install MCNV2 (R)
--------------------

.. code-block:: r

   if (!require("devtools")) install.packages("devtools")

   devtools::install_github(
     "JacquemontLab/MCNV2-Mendelian-CNV-Validation",
     ref = "package-integration",
     dependencies = TRUE
   )

**Alternative:** Install from local tarball

.. code-block:: r

   install.packages("/path/to/MCNV2_0.1.0.tar.gz", repos = NULL, type = "source")

2) Python requirements
----------------------

MCNV2 uses a Python virtual environment to ensure reproducibility.

**Create and configure the environment:**

.. code-block:: r

   library(MCNV2)
   # Creates 'r-MCNV2' virtualenv and installs dependencies (polars, etc.)
   MCNV2::setup_python_env(envname = "r-MCNV2")

Dependencies are defined in: ``inst/python/requirements.txt``

**Activate the environment:**

.. code-block:: r

   library(reticulate)
   use_virtualenv("r-MCNV2", required = TRUE)

**Verify configuration:**

.. code-block:: r

   py_config()

**Optional: check that key packages are available:**

.. code-block:: r

   py_run_string("import polars; print(polars.__version__)")

3) Install bedtools
-------------------

bedtools is required for genomic overlap calculations.

**macOS (Homebrew):**

.. code-block:: bash

   brew install bedtools

**Ubuntu/Debian:**

.. code-block:: bash

   sudo apt-get install bedtools

**Conda (any platform):**

.. code-block:: bash

   conda install -c bioconda bedtools

**Verify installation:**

.. code-block:: bash

   bedtools --version

4) Launch MCNV2
---------------

After completing the installation:

.. code-block:: r

   library(reticulate)
   use_virtualenv("r-MCNV2", required = TRUE)

   library(MCNV2)
   MCNV2::launch(
     bedtools_path = Sys.which("bedtools"),
     results_dir = "~/mcnv2_results"
   )

This will open the interactive Shiny application for CNV validation and annotation.

.. tip::

   For batch processing and reproducible pipelines, see the :doc:`../tutorials/cli_tutorial`.
