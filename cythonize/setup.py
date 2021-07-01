from setuptools import Extension, setup
from Cython.Build import cythonize

setup(
    name="vaspwfc cythonized",
    ext_modules=cythonize(Extension(
	"vaspwfc", sources=["VaspBandUnfolding/cythonize/vaspwfc.pyx"], language="c++", build_dir="VaspBandUnfolding/cythonize", extra_compile_args=["-std=c++11"], extra_link_args=["-std=c++11"]), annotate=True)
)
