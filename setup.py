from setuptools import setup, find_packages
import thermopy


setup(name="thermopy",
      version=thermopy.__version__,
      description=str('Some utilities for Thermodynamics and Thermochemistry'
                      'using the NASA 9 polynomials'),
      author="Felipe M. Vieira",
      author_email="fmv1992@gmail.com",
      url="",
      license="GPL",
      packages=find_packages(),
      include_package_data=True,
      package_data={
          'database': ['burcat_thr.xml', 'nasa9polynomials.xml'],
          'documentation': ['thermopy050_documentation.pdf',
                            'thermopy050_overview.pdf']},
      install_requires=['scipy>=0.6.0', 'numpy>=1.2.1'],
      test_suite='nose.collector',
      tests_require=['nose'],
      zip_safe=False,
      keywords='thermodynamics, properties estimation',
      classifiers=[
          'Development Status :: Alpha',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: GPL',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3.4',
          'Topic :: Scientific/Engineering :: Chemistry',
          'Topic :: Scientific/Engineering :: Physics',
          'Topic :: Software Development :: Libraries :: Python Modules'
      ]
      )
