from setuptools import setup
from Cython.Build import cythonize

setup(
    name="vaspwfc cythonized",
    ext_modules=cythonize("VaspBandUnfolding/cythonize/vaspwfc.pyx"),
    zip_safe=False,
)
