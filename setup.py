#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup
import pycoQC as pqc

setup(
    name = pqc.__name__,
    version = pqc.__version__,
    description = pqc.__description__,
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
        'scipy==1.3.0',
        'pandas==0.24.2',
        'plotly==3.10.0',
        'jinja2==2.10.1',
        'h5py==2.9.0',
        'tqdm==4.32.1'],
    packages = [pqc.__name__],
    package_dir = {pqc.__name__: pqc.__name__},
    package_data = {pqc.__name__: ['templates/*']},
    entry_points = {
        'console_scripts': [
            'pycoQC=pycoQC.__main__:main_pycoQC',
            'Fast5_to_seq_summary=pycoQC.__main__:main_Fast5_to_seq_summary']}
)
