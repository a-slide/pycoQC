#!python3
# -*- coding: utf-8 -*-

from setuptools import setup

# Define package info
with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
    name = 'pycoQC',
    description = 'PycoQC computes metrics and generates interactive QC plots for Oxford Nanopore technologies sequencing data',
    version = '2.5.2',
    long_description = long_description,
    long_description_content_type='text/markdown',
    url = 'https://github.com/a-slide/pycoQC',
    author = 'Adrien Leger & Tommaso Leonardi',
    author_email = 'aleg@ebi.ac.uk',
    license = 'GPLv3',
    python_requires ='>=3.6',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3'],
    install_requires = [
        'numpy>=1.19',
        'scipy>=1.5',
        'pandas>=1.1',
        'plotly==4.1.0',
        'jinja2>=2.10',
        'h5py>=3.1',
        'tqdm>=4.54',
        'pysam>=0.16'],
    packages = ['pycoQC'],
    package_dir = {'pycoQC': 'pycoQC'},
    package_data = {'pycoQC': ['templates/*']},
    entry_points = {
        'console_scripts': [
            'pycoQC=pycoQC.__main__:main_pycoQC',
            'Fast5_to_seq_summary=pycoQC.__main__:main_Fast5_to_seq_summary',
            'Barcode_split=pycoQC.__main__:main_Barcode_split']}
)
