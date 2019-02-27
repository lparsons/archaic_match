#!/usr/bin/env python
# -*- coding: utf-8 -*-


from setuptools import setup

version = {}
with open("archaic_match/version.py") as fp:
    exec(fp.read(), version)

setup(
    name='archaic_match',
    version=version['__version__'],
    packages=['archaic_match'],
    entry_points={
        'console_scripts': [
            'archaic_match=archaic_match.__main__:main',
            'column_replace=archaic_match.column_replace:main'
        ]
    },
    install_requires=[
        'pandas',
        'scikit-allel>=1.2.0'
    ]
)
