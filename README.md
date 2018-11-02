# pycoQC-2

[![DOI](https://zenodo.org/badge/94531811.svg)](https://zenodo.org/badge/latestdoi/94531811)
[![GitHub license](https://img.shields.io/github/license/a-slide/pycoQC.svg)](https://github.com/a-slide/pycoQC/blob/master/LICENSE)
[![Language](./pictures/language-Python3-brightgreen.svg)](https://www.python.org/)

---

**PycoQC-2 is a Pure Python 3 package for Jupyter Notebook, computing metrics and generating simple QC plots from the sequencing summary report generated by Oxford Nanopore technologies Albacore basecaller**

As opposed to more exhaustive QC programs for nanopore data, pycoQC is very fast as it relies entirely on the *sequencing_summary.txt* file generated by ONT Albacore basecaller. Consequently, pycoQC will only provide metrics at read level metrics (and not at base level). The package supports 1D and 1D2 runs analysed with Albacore.

PycoQC requires the following fields in the sequencing.summary.txt file:

* 1D run => **read_id**, **run_id**, **channel**, **start_time**, **sequence_length_template**, **mean_qscore_template**
* 1D2 run =>**read_id**, **run_id**, **channel**, **start_time**, **sequence_length_2d**, **mean_qscore_2d**

In addition it will try to get the following optional fields if they are available:

* **calibration_strand_genome_template**, **barcode_arrangement**

# Gallery

*Click on picture to see online interactive version.* 

[![reads_len_1D_example](./pictures/reads_len_1D_example.png)](https://plot.ly/~aleg/2/distribution-of-read-length/)

[![](./pictures/summary.png)](https://plot.ly/~aleg/8/)

[![reads_qual_len_2D_example](./pictures/reads_qual_len_2D_example.png)](https://plot.ly/~aleg/3/mean-read-quality-per-sequence-length/)

[![channels_activity](./pictures/channels_activity.png)](https://plot.ly/~aleg/4/output-per-channel-over-experiment-time/)

[![output_over_time](./pictures/output_over_time.png)](https://plot.ly/~aleg/5/output-over-experiment-time/)

[![barcode_counts](./pictures/barcode_counts.png)](https://plot.ly/~aleg/7/percentage-of-reads-per-barcode/)

[![qual_over_time](./pictures/qual_over_time.png)](https://plot.ly/~aleg/6/mean-read-quality-over-experiment-time/)

# Installation

Ideally, before installation, create a clean **Python 3** virtual environment to deploy the package. **Python 2 is not supported**

For example you can use virtualenvwrapper (see http://www.simononsoftware.com/virtualenv-tutorial-part-2/).

## Dependencies

pycoQC relies on a few robustly maintained third party libraries (numpy, scipy, plotly, pandas). The correct versions of the packages are installed together with the software when using pip.

## Option 1: Direct installation with pip from github (recommended)

Install the package with pip3. Python dependencies are automatically installed.

`pip3 install git+https://github.com/a-slide/pycoQC.git`

To update the package:

`pip3 install git+https://github.com/a-slide/pycoQC.git --upgrade`


## Option 2: Clone the repository and install locally in develop mode

With this option, the package will be locally installed in “editable” or “develop” mode. This allows the package to be both installed and editable in project form. This is the recommended option if you wish to participate to the development of the package. Python dependencies will be automatically installed.

`git clone https://github.com/a-slide/pycoQC.git`

`cd pycoQC`

`chmod u+x setup.py`

`pip3 install -e ./`

# Usage in Jupyter Notebook

## Demo notebook

An online live usage notebook served by MyBinder is available to familiarize with the package API:

[![mybinder](./pictures/launch-mybinder-red.svg)](https://mybinder.org/v2/gh/a-slide/pycoQC/dev?filepath=tests%2FpycoQC_usage.ipynb) (Can take a few minutes)

A static html version of the same notebook is also available if you experience any issue with the live version:

[![html](./pictures/static-html-blue.svg)](http://htmlpreview.github.io/?https://github.com/a-slide/pycoQC/blob/dev/tests/pycoQC_usage.html)

## Running your own notebook locally or remotely

If you want to run pycoQC interactively in Jupyter you need to install Jupyter manually.

If you installed pycoQC in a virtual environment, be carefull to install jupyter notebook in the same virtual environment.

`pip3 install notebook`

Launch the notebook in a shell terminal

`jupyter notebook`

If it does not autolaunch your web browser, open manually the following URL http://localhost:8888/tree

From Jupyter home page you can navigate to the directory you want to work in. Then, create a new Python3 Notebook.


# Shell standalone interface

pycoQC 2 now has the ability to generate an HTML report containing all the plots.

TODO


# Note to power-users and developers

Please be aware that pycoQC is an experimental package that is still under development. It was tested under Linux Ubuntu 16.04 and in an HPC environment running under Red Hat Enterprise 7.1.

You are welcome to contribute by requesting additional functionalities, reporting bugs or by forking and submitting pull requests

Thank you

### Authors

* Adrien Leger - aleg {at} ebi.ac.uk
