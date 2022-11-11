from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy

setup(
    name="vaspwfc cythonized",
    ext_modules=cythonize(Extension(
	"vaspwfc", sources=["VaspBandUnfolding/cythonize/vaspwfc.pyx"], language="c++", build_dir="VaspBandUnfolding/cythonize", extra_compile_args=["-std=c++11"], extra_link_args=["-std=c++11"], include_dirs=[numpy.get_include()]), annotate=True, compiler_directives={'language_level': "3"})
)
