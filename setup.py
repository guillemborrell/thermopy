from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages

setup(name="thermopy",
      version="0.2",
      description="Some utilities for Thermodynamics and Thermochemistry",
      author="Guillem Borrell i Nogueras",
      author_email="guillem@torroja.dmt.upm.es",
      url="http://torroja.dmt.upm.es/guillem/blog/",
      packages = find_packages(),
      include_package_data = True,
      install_requires = ['scipy>=0.6.0','numpy>=1.2.1'],
      zip_safe = False
      )
