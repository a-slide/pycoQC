#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""       
  ___              ___   ___ 
 | _ \_  _ __ ___ / _ \ / __|
 |  _/ || / _/ _ \ (_) | (__ 
 |_|  \_, \__\___/\__\_\\___|
      |__/      
                                __   __     ___ 
 /\  _| _. _ _   |   _ _  _ _    _) /  \ /|   / 
/--\(_|| |(-| )  |__(-(_)(-|    /__ \__/  |  /  
                      _/                        
"""

from setuptools import setup

setup(
    name='pycoQC',
    version='1.0.dev3',
    description='Python 3 package for Jupyter Notebook, computing metrics and generating simple QC plots from Oxford Nanopore technologies (ONT) Albacore basecaller',
    long_description=open('README.rst', 'r').read(),
    url='https://github.com/a-slide/pycoQC',
    author='Adrien Leger',
    author_email='aleg@ebi.ac.uk',
    license='GPLv3',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',],
    install_requires=[
        'numpy>=1.13.0',
        'pandas>=0.20.0',
        'matplotlib>=2.0.0',
        'seaborn>= 0.7.0',
        'notebook>=4.0.0',],
    packages=["pycoQC"],
    package_dir={'pycoQC': 'pycoQC'},
    package_data={'pycoQC': ['data/sequencing_summary.txt', "test_pycoQC.ipynb"]},

)
