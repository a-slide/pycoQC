#!bash
# -*- coding: utf-8 -*-

set -e

echo "print dir content"
ls
ls dist

echo "Build noarch package..."
conda build meta.yaml --python 3.6 --numpy 1.1 --output-folder conda_build --debug

echo "Deploying to Anaconda.org..."
anaconda -t $ANACONDA_TOKEN upload conda_build/**/pycoqc-*.tar.bz2

echo "Successfully deployed to Anaconda.org."
exit 0
