# thermopy
Python library for thermodynamics and other handy tools.


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

Inside thermopy directory execude:
```
python3 thermopy/setup.py test
```

## Installing

Once the tests ran fine install it:
```
python3 thermopy/setup.py install
```


-----BEGIN PGP SIGNED MESSAGE-----

Hash: SHA256

sha1sum of the files:
02cec2818e3c9c94773d9aa05e006239bbd8358e  ./authors
bf21a9e8fbc5a3846fb05b4fa0859e0917b2202f  ./.cache/v/cache/lastfailed
f2f7546678f0bc91f0c697b450cb7a409d73b522  ./license
6d659e8b9f441cb3f7bf87778925b3a263313964  ./documentation/thermopy050_documentation.pdf
7a168489b6a7a4330be4d35eaeae69b815fe8988  ./documentation/thermopy050_overview.pdf
b73c8eed39cb2ef27ebc939804059039c1dfb0b9  ./tests/__init__.py
7c935759c9f198665b45773fbb27badfbf4f9b95  ./tests/test_burcat.py
f84a2e1703d04cdbeb2465984397b295ee76cd42  ./tests/test_iapws.py
d54fd80fc2ab49c2c7d7dd725a1c866b4f0783e8  ./tests/test_units.py
ae3713666631c055eeb9dabbdea89531ab408a9b  ./thermopy/burcat.py
b08716195c4f58e90d2a736f322a63ba91656fef  ./thermopy/units.py
f68baaa715b111a7e6cd1cb7f51e24c28de23cb3  ./thermopy/constants.py
b73c8eed39cb2ef27ebc939804059039c1dfb0b9  ./thermopy/__init__.py
e7a0d270f3d0ac91392f28599181584f1481ea31  ./thermopy/nasa9polynomials.py
1e5152232b0d8a390bce9e5e117c6d6fd5ffcc37  ./thermopy/iapws.py
92e35e5bb7d55016152127902e8401207da91e37  ./__init__.py
e3259d616a62ed090b4fc6e384ff5fdb46154210  ./databases/nasa9polynomials.xml
9ec73c56fbd30560ab73c0028d477dfa441b35ac  ./databases/burcat_thr.xml
e68f711a958575dac5aadf5e47f8de53e80e6621  ./setup.py

-----BEGIN PGP SIGNATURE-----

Version: GnuPG v1


iQEcBAEBCAAGBQJWPqRzAAoJED+RfayjyZfdUQwIAMS1/qEJ5zMqGxDcuvtqF9Dn
N/I5Q1LDByOmTyfiT9F52LQ4/ozUB0IAb1qsJEMQ4tA0F7sUMAKeN2Jf9oCiZJN6
WF4D4Xp9nruR/Mrm9cAodjkDjiVk+W/PslMVys5UvxYmj3dxFkgCb/8XNv7FyaSq
Q+UQjriDTSS4Ip/DKmnHlGKR6KHUtCh6XHnMmcpUj18S0Rqbfud2VaK5+m7fw0QC
ClLoQ7A+u1bCjN3fl85YGXzVbFlRBGQieS4irDzaze1bXgpTQZcVfsA0vuVkN1W7
agUoLdDF911/jcjo72OsuKtTLGfTF/TF4aSwaro62VjHNl+DZxNC/8EwSJYxf9g=
=seJo

-----END PGP SIGNATURE-----


key: https://pgp.mit.edu/pks/lookup?op=vindex&search=0x7BCA19BB0E69E45D
