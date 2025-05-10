# setup.py
from setuptools import setup
from Cython.Distutils import build_ext


setup(
    name="ABPy",
    packages=["ABPy"],
    cmdclass={'build_ext': build_ext},
    include_dirs=['./src/']
)