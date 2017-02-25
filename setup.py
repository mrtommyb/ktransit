#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from numpy.distutils.core import setup, Extension



if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()


# Define the Fortran extension.
tmodtom = Extension("ktransit._tmodtom", ["ktransit/tmodtom.f"])

setup(
    name="ktransit",
    description='A simple exoplanet transit modeling tool in python.',
    url="https://github.com/mrtommyb/ktransit",
    version="0.2.5",
    author="Tom Barclay",
    author_email="tom@tombarclay.com",
    packages=["ktransit", ],
    package_dir={'ktransit':'ktransit'},
    ext_modules=[tmodtom, ],
    classifiers=[
        # "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
