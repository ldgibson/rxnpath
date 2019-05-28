from __future__ import absolute_import, division, print_function
from os.path import join as pjoin

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
_version_major = 0
_version_minor = 1
_version_micro = ''  # use '' for first of series, number for 1 and above
_version_extra = 'dev2'
_version_extra = ''  # Uncomment this for full releases

# Construct full version string from these.
_ver = [_version_major, _version_minor]
if _version_micro:
    _ver.append(_version_micro)
if _version_extra:
    _ver.append(_version_extra)

__version__ = '.'.join(map(str, _ver))

CLASSIFIERS = ["Development Status :: 3 - Alpha",
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: MIT License",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering"]

# Description should be a one-liner:
description = "rxnpath: a simple python package for generating reaction diagrams."
# Long description will go up on the pypi page
long_description = """

rxnpath
========
rxnpath is a python package for generating custom reaction diagrams.

License
=======
``rxnpath`` is licensed under the terms of the MIT license. See the file
"LICENSE" for information on the history of this software, terms & conditions
for usage, and a DISCLAIMER OF ALL WARRANTIES.

All trademarks referenced herein are property of their respective holders.

Copyright (c) 2019--, Luke D Gibson, The University of Washington
Department of Chemical Engineering.
"""

NAME = "rxnpath"
MAINTAINER = "Luke D Gibson"
MAINTAINER_EMAIL = "ldgibson@uw.edu"
DESCRIPTION = description
LONG_DESCRIPTION = long_description
URL = "http://github.com/ldgibson/rxnpath"
DOWNLOAD_URL = ""
LICENSE = "MIT"
AUTHOR = "Luke D Gibson"
AUTHOR_EMAIL = "ldgibson@uw.edu"
PLATFORMS = "OS Independent"
MAJOR = _version_major
MINOR = _version_minor
MICRO = _version_micro
VERSION = __version__
PACKAGE_DATA = {'rxnpath': [pjoin('data', '*')]}
REQUIRES = ['adjustText>=0.7.3', 'matplotlib>=3.1.0', 'networkx>=2.3', 'numpy>=1.16.3']
