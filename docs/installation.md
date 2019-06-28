# Installation

## Create a clean virtual environment

Ideally, before installation, create a clean **python3.5+** virtual environment to deploy the package. **Python 2 is not supported**.
For example one can use conda or virtualenvwrapper.

With [virtualenvwrapper](https://virtualenvwrapper.readthedocs.io/en/latest/install.html):

```bash
mkvirtualenv pycoQC -p python3.6
```

With [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

```bash
conda create -n pycoQC python=3.6
```

## Dependencies

pycoQC relies on a few robustly maintained third party libraries listed below. The correct versions of the packages are installed together with the software when using pip.

* numpy>=1.13
* scipy>=1.1
* pandas>=0.23
* plotly>=3.4
* jinja2>=2.10
* h5py>=2.8.0
* tqdm>=4.23'

## Option 1: Installation with pip from pypi

```bash
# First installation
pip install pycoQC

# Update to last version
pip install pycoQC --upgrade
```

## Option 2: Installation with conda from Anacounda cloud

In a previously created conda environment

```bash
# First installation
conda install -c aleg pycoqc

# Update to last version
conda update -c aleg pycoqc
```

## Option 3: Installation with pip from Github

To get the last stable (master) or bleeding edge (dev) version

```bash
# First installation
pip install git+https://github.com/a-slide/pycoQC.git

# First installation bleeding edge
pip install git+https://github.com/a-slide/pycoQC.git@dev

# Update to last version
pip install git+https://github.com/a-slide/pycoQC.git --upgrade
```

## Option 4: Clone the repository and install locally in develop mode

With this option, the package will be locally installed in *editable* or *develop mode*. This allows the package to be both installed and editable in project form. This is the recommended option if you wish to modify the code and/or participate to the development of the package (see [contribution guidelines](contributing.md)).

```bash
# Clone stable repo locally
git clone https://github.com/a-slide/pycoQC.git

# bleeding edge branch
git clone --branch dev https://github.com/a-slide/pycoQC.git

# Enter in repo directory
cd pycoQC

# Make setup.py executable
chmod u+x setup.py

# Install with pip3
pip3 install -e ./
```
