# setup.py
from setuptools import find_packages, setup, Extension
from Cython.Distutils import build_ext

ext_modules=[
    Extension(
        "ABPy.triangular_lattice", # where the .so will go
        ["./bindings/triangular_lattice.pyx"]
    )
]

setup(
    name="ABPy",
    packages=["ABPy"],
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules,
    include_dirs=['./src/']
)