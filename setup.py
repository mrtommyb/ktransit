#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from numpy.distutils.core import setup, Extension


if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()


# First, make sure that the f2py interfaces exist.
interface_exists = os.path.exists("ktransit/tmodtom.pyf")
if "interface" in sys.argv or not interface_exists:
    # Generate the Fortran signature/interface.
    cmd = "cd k-transit;"
    cmd += "f2py tmodtom.f -m _tmodtom -h tmodtom.pyf"
    cmd += " --overwrite-signature"
    os.system(cmd)
    if "interface" in sys.argv:
        sys.exit(0)

# Define the Fortran extension.
tmodtom = Extension("ktransit._tmodtom", ["ktransit/tmodtom.pyf", "ktransit/tmodtom.f"])

setup(
    name="ktransit",
    url="https://github.com/mrtommyb/ktransit",
    version="0.2.2",
    author="Tom Barclay",
    author_email="tom@tombarclay.com",
    description="",
    long_description="",
    packages=["ktransit", ],
    ext_modules=[tmodtom, ],
    classifiers=[
        # "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
