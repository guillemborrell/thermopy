# thermopy3
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

## Testing

Inside thermopy directory execute:
```
python3 thermopy/setup.py test
```

## Installing

Once the tests ran fine install it:
```
python3 thermopy/setup.py install
```

## Changelog:

v0.5.2:

    - Changed the names from thermopy to thermopy3 because the former was already in use in pypi.
    
    - IAPWS is now on SI and (molar basis) instead of its native kJ, kg units.

v0.5.1:

	- Fix error with relative import of xml databases

v0.5.0:

	- First release from v0.4.0
	
	- Ported from python2.7 to python3.4
