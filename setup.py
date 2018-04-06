#!/usr/bin/env python
# -*- coding: utf-8 -*-


from setuptools import setup

setup(
    name='archaic_match',
    version='0.1.0',
    packages=['archaic_match'],
    entry_points={
        'console_scripts': [
            'archaic_match=archaic_match.__main__:main'
        ]
    },
    install_requires=[
        'sqlalchemy',
        'numpy',
        'scipy',
        'matplotlib',
        'seaborn',
        'pandas',
        'scikit-learn',
        'h5py!=2.7.0,!=2.7.1',  # These versions issue warnings when imported
        'numexpr',
        'bcolz',
        'zarr',
        'dask',
        'cytoolz',
        'scikit-allel'
    ]
        )
