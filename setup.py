from setuptools import setup, find_packages
import thermopy

my_long_description = str(
"""
Python library for thermodynamics and other handy tools.


Thermodynamics (all these properties as function of temperature for thousands
of compounds):

	- Specific heat capacity

	- Enthalpy

	- Entropy


Temperature independent data:

	- Molecular weight

	- Enthalpy of formation


and much more.

For water pressure is also an input (higher accuracy).

Modelling of chemical reactions is also present. Main features:

	- Equilibrium constant as a function of temperature

	- Heat of reaction as a function of temperature


Handy tools:

	- Units conversion module

	- Hundreds of physical constants


See the documentation for further details and examples.
""")


setup(name="thermopy",
      version=thermopy.__version__,
      description='Python package for thermodynamic calculations and units '
                  'conversion.',
      long_description = my_long_description,
      author="Felipe M. Vieira",
      author_email="fmv1992@gmail.com",
      url="github: https://github.com/guillemborrell/thermopy",
      license="GPL",
      packages=find_packages(),
      include_package_data=True,
      data_files=[('databases', ['databases/burcat_thr.xml',
                                 'databases/nasa9polynomials.xml'])],
      install_requires=['scipy>=0.6.0', 'numpy>=1.2.1'],
      test_suite='nose.collector',
      tests_require=['nose'],
      zip_safe=False,
      keywords='thermodynamics, properties estimation',
      # full listing on https://pypi.python.org/pypi?%3Aaction=list_classifiers
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 or later'
          ' (GPLv3+)',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3.4',
          'Topic :: Scientific/Engineering :: Chemistry',
          'Topic :: Scientific/Engineering :: Physics',
          'Topic :: Software Development :: Libraries :: Python Modules'
      ]
      )
