# -*- coding: utf-8 -*-

# Define self package variable
__version__ = '2.2.2'
__all__ = ["pycoQC", "common"]
__description__="""
PycoQC computes metrics and generates interactive QC plots for Oxford Nanopore technologies sequencing data
"""
# Collect info in a dictionnary for setup.py
setup_dict = {
    "name": __name__,
    "version": __version__,
    "description": __description__,
    "url": "https://github.com/a-slide/pycoQC",
    "author": 'Adrien Leger / Tommaso Leonardi',
    "author_email": 'aleg {at} ebi.ac.uk / tom {at} tleo.io',
    "license": 'GPLv3',
    "python_requires":'>=3.3',
    "classifiers": [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3'],
    "install_requires": [
        'numpy>=1.13',
        'scipy>=1.1',
        'pandas>=0.23',
        'plotly>=3.4',
        'jinja2>=2.10',
        'h5py>=2.8.0',
        'tqdm>=4.23'],
    "packages": [__name__],
    "package_dir": {__name__: __name__},
    "package_data": {__name__: ['templates/*']},
    "entry_points": {
        'console_scripts': [
            'pycoQC=pycoQC.cli:main_pycoQC',
            'Fast5_to_seq_summary=pycoQC.cli:main_Fast5_to_seq_summary'],
    }
}
