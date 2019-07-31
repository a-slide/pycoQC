#!python3
# -*- coding: utf-8 -*-

from setuptools import setup

# Define package info
name = "pycoQC"
version = "2.3.1.6"
description = "PycoQC computes metrics and generates interactive QC plots for Oxford Nanopore technologies sequencing data"
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name = name,
    description = description,
    version = version,
    long_description = long_description,
    long_description_content_type="text/markdown",
    url = "https://github.com/a-slide/pycoQC",
    author = 'Adrien Leger & Tommaso Leonardi',
    author_email = 'aleg@ebi.ac.uk',
    license = 'GPLv3',
    python_requires ='>=3.5',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3'],
    install_requires = [
        'numpy==1.16.4',
        'scipy==1.2.0',
        'pandas==0.24.2',
        'plotly==3.10.0',
        'jinja2==2.10.1',
        'h5py==2.9.0',
        'tqdm==4.32.1'],
    packages = [name],
    package_dir = {name: name},
    package_data = {name: ['templates/*']},
    entry_points = {
        'console_scripts': [
            'pycoQC=pycoQC.__main__:main_pycoQC',
            'Fast5_to_seq_summary=pycoQC.__main__:main_Fast5_to_seq_summary']}
)
