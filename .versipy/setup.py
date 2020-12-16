#!python3
# -*- coding: utf-8 -*-

from setuptools import setup

# Define package info
with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
    name = '__package_name__',
    description = '__package_description__',
    version = '__package_version__',
    long_description = long_description,
    long_description_content_type='text/markdown',
    url = '__package_url__',
    author = '__author_name__',
    author_email = '__author_email__',
    license = '__package_licence__',
    python_requires ='__minimal_python__',
    classifiers = [
        '__classifiers_1__',
        '__classifiers_2__',
        '__classifiers_3__',
        '__classifiers_4__',
        '__classifiers_5__'],
    install_requires = [
        '__dependency_1__',
        '__dependency_2__',
        '__dependency_3__',
        '__dependency_4__',
        '__dependency_5__',
        '__dependency_6__',
        '__dependency_7__',
        '__dependency_8__'],
    packages = ['__package_name__'],
    package_dir = {'__package_name__': '__package_name__'},
    package_data = {'__package_name__': ['templates/*']},
    entry_points = {
        'console_scripts': [
            '__entry_point_1__',
            '__entry_point_2__',
            '__entry_point_3__']}
)
