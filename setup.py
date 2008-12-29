from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages

setup(name="thermopy",
      version="0.1",
      description="Some utilities for Thermodynamics and Thermochemistry",
      author="Guillem Borrell i Nogueras",
      author_email="guillem@torroja.dmt.upm.es",
      url="http://torroja.dmt.upm.es/guillem/blog/",
      packages = find_packages(),
      package_data = {
        '': ['*.xml', '*.rst'],
        },
      install_requires = ['numpy>=1.2.1','scipy>=0.6.0'],
      )
