# run
#      python setup.py build_ext --inplace

from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
setup(
    ext_modules = cythonize(Extension('croutines',["croutines.pyx"]), language_level=3),
    include_dirs=[numpy.get_include()]
)
