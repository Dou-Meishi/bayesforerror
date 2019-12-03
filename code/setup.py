#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from distutils.core import setup, Extension
from Cython.Build import cythonize

setup(
        ext_modules = cythonize(Extension(
            'c_priors',
            sources = ['c_priors.pyx'],
            include_dirs = [np.get_include()],
            # extra_compile_args = ['-ffast-math']
        ),  annotate=True)
)
