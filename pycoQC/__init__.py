# -*- coding: utf-8 -*-

# Define self package variable
__version__ = '1.0a1'
__all__ = ["pycoQC", "pycoQC_fun"]
__doc__='Python 3 package for Jupyter Notebook, computing metrics and generating simple QC plots from Oxford Nanopore technologies (ONT) Albacore basecaller'
__author__= 'Adrien Leger'
__email__ = 'aleg@ebi.ac.uk'
__url__ = "https://github.com/a-slide/pycoQC"
__licence__ = 'GPLv3'
__classifiers__ = [
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',]
__install_requires__ = ['numpy>=1.13.0', 'pandas>=0.20.0', 'matplotlib>=2.0.0', 'seaborn>= 0.7.0', 'notebook>=4.0.0']
__package_data__ =  ['data/sequencing_summary.txt', 'data/sequencing_1dsq_summary.txt', "test_pycoQC.ipynb"]


# Collect info in a dictionnary for setup.py
from collections import OrderedDict
setup_dict = OrderedDict()

setup_dict["name"]=__name__
setup_dict["version"]=__version__
setup_dict["description"]=__doc__
setup_dict["url"]=__url__
setup_dict["author"]= __author__
setup_dict["author_email"]=__email__
setup_dict["license"]=__licence__
setup_dict["classifiers"]= __classifiers__
setup_dict["install_requires"]= __install_requires__
setup_dict["packages"]=[__name__]
setup_dict["package_dir"]={__name__: __name__}
setup_dict["package_data"]={__name__: __package_data__}
