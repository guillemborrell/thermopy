## Note:

This repo just received https://github.com/fmv1992/thermopy3. There is an
ongoing process to join the old thermopy (python2) with thermopy (python3).
Soon there will be only thermopy supporting python3.

# thermopy / thermopy3 (see note above)

Python library for thermodynamics and other handy tools. The library was not extensively tested
with python2.

Thermodynamics (all these properties as function of temperature for thousands of compounds):

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

## Installing
Make sure you have both `numpy` and `scipy` installed.

Then install it:
```
python3 thermopy/setup.py install
```
## Testing

Inside thermopy directory execute:
```
python3 thermopy/setup.py test
```

## Changelog:

v0.5.4:

    - Merge complete.

    - Documentation now using sphinx instead of pdf files.

    - Code compliance and readability improved.

    - Code and logic uniformized.

v0.5.3 (ongoing):
	
	- Merged with original thermopy and discarded the 'thermopy3' name.

	- Added meaningful docstrings to every package, module, class, method and function.

	- Uniformized docstrings to comply with google docstrings style.

	- Improved compliance with PEP 257 and PEP 8.

	- Migrated documentation to sphinx (pdf will no longer be available).

v0.5.2:

    - Changed the names from thermopy to thermopy3 because the former was already in use in pypi.
    
    - IAPWS is now on SI and (molar basis) instead of its native kJ, kg units.

v0.5.1:

	- Fix error with relative import of xml databases.

v0.5.0:

	- First release from v0.4.0.
	
	- Ported from python2.7 to python3.4.

## TODO

	- Transform IAPWS module into molar basis to be consistent with burcat and nasa9polynomials.

	- Fix the setting of a compound when it is ambiguous:
		>>> csbr = db.set_compound('LYQFWZFBNBDLEO-UHFFFAOYSA-M')
		Traceback (most recent call last):
		  File "<stdin>", line 1, in <module>
		  File "/home/monteiro/cloud/cloud_work/fmv1992_github/thermopy/thermopy/nasa9polynomials.py", line 501, in set_compound
			('Ga2O', 'gallium;oxygen(2-)')
		Exception: ("The compound 'LYQFWZFBNBDLEO-UHFFFAOYSA-M' you are trying to set is not unique: CsBr", 'CsBr(cr)')

	- Put examples again in docstrings; need to update.

	- Uniformize code (imports and the like) and improve compliance with PEP8.

	- Update docstrings to conform to PEP 257.

	- Update documentation to use the package's name after merge: 'thermopy'.

	- Do more testing. It looks like in some cases there are some import errors
	  (see 'EXAMPLE 02: Hydrazine ’messing around’ example.' in documentation). 

	- Make the PDFs 'copy and pastable'.

	- Increase the testing coverage.

    - Implement units for every output so output is not a number but a dimension (e.g. 1 J/kg instead of 1).

