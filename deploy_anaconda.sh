#!bash
# -*- coding: utf-8 -*-

set -e

echo "Set up conda package manager"
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh --quiet
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no --set anaconda_upload no
conda update -q conda

echo "Install packages needed for package build and upload"
conda install -q python=3.6 conda-build anaconda-client ripgrep conda-verify

echo "compile package from setup.py"
python setup.py sdist

echo "Build noarch package..."
conda build meta.yaml --python 3.6 --numpy 1.1 --output-folder conda_build  -c bioconda -c conda-forge

echo "Deploying to Anaconda.org..."
anaconda -v -t $1 upload conda_build/**/*.tar.bz2

echo "Successfully deployed to Anaconda.org."
exit 0
