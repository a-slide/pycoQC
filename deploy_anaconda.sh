#!bash

set -e

echo "Build recipe from pypi package..."
conda skeleton pypi pycoQC --output-dir conda_build --noarch-python

echo "Build noarch package..."
conda build conda_build/pycoqc/meta.yaml --output-folder conda_build

echo "Deploying to Anaconda.org..."
anaconda -t $ANACONDA_TOKEN upload conda_build/**/pycoqc-*.tar.bz2

echo "Successfully deployed to Anaconda.org."
exit 0
