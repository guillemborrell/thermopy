# thermopy3
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
<<<<<<< HEAD
```

## Cryptographic signatures:

### v0.5.1:

#### github

-----BEGIN PGP SIGNED MESSAGE-----

Hash: SHA256

e581449592c0e6d410b176688b1c3396b0abcc6f  ./documentation/thermopy051_documentation.pdf

88c1aa5ecf7a44a666962c243c33e9c3d4581c46  ./documentation/thermopy051_overview.pdf

d881699e8ab6ee680f6d243d6736bba88a252dca  ./thermopy/burcat.py

b08716195c4f58e90d2a736f322a63ba91656fef  ./thermopy/units.py

f68baaa715b111a7e6cd1cb7f51e24c28de23cb3  ./thermopy/constants.py

e555eefb44497199719d371afb70538c5ecf2f34  ./thermopy/\_\_init\_\_.py

cf04cd22042b126795b8e4652f7cd21732cd863a  ./thermopy/nasa9polynomials.py

1e5152232b0d8a390bce9e5e117c6d6fd5ffcc37  ./thermopy/iapws.py

1f577957e736579b8a58bede5824bb5a673b9ec6  ./\_\_init\_\_.py

f2f7546678f0bc91f0c697b450cb7a409d73b522  ./license.txt

02cec2818e3c9c94773d9aa05e006239bbd8358e  ./authors.txt

7f43802016fdbdd86c4fa381e2c33525fd04d140  ./MANIFEST.in

e3259d616a62ed090b4fc6e384ff5fdb46154210  ./databases/nasa9polynomials.xml

9ec73c56fbd30560ab73c0028d477dfa441b35ac  ./databases/burcat_thr.xml

c6bbf3fd7fac4108000b10367fa05ed0a247e35b  ./setup.py

b73c8eed39cb2ef27ebc939804059039c1dfb0b9  ./test/\_\_init\_\_.py

7c935759c9f198665b45773fbb27badfbf4f9b95  ./test/test_burcat.py

f84a2e1703d04cdbeb2465984397b295ee76cd42  ./test/test_iapws.py

d54fd80fc2ab49c2c7d7dd725a1c866b4f0783e8  ./test/test_units.py

37af826fab2067fcd4044ac4c17093f88c1d5e64  ./README.md

-----BEGIN PGP SIGNATURE-----

Version: GnuPG v1

iQEcBAEBCAAGBQJWP/tkAAoJED+RfayjyZfdVu0H/A2dt2uKurN03rqiiZSIa4fa
rYazGO4Fz7MbAZS5ybRJb+55TXnfbCFlrT2bZrY9FWUkN9QGTtAHH587PwIvpefs
SUQ2rF04VcRnXNvuN6PmGqRlhtC9iFXa5yTFHkGhuFlU0hR15OxXg6pxrKR+a68V
kco8rnX34wMQwxF06BCW9VlNnV3+EqYPqv78b2RuYt4wAdo17JhY25WyeK6jIjew
5txjOggr+mTkVg7DNRORgibmMmHXY6S0GpFhr/qZACj1R/x1yo3PslofQYI5LdPC
A18AAw/wd3R+F0F6xVCBuv8pbUAfyRZ7GkcvKdYPCT49gYF42rLwJUAS8TrINMA=
=QIY6

-----END PGP SIGNATURE-----

#### pypi:

sha1sum:	0860fea268e8e89fc6e9b697b78236a32208ce63  thermopy3-0.5.1.tar.gz

-----BEGIN PGP SIGNATURE-----

Version: GnuPG v1

iQEcBAABCAAGBQJWP/uaAAoJED+RfayjyZfdcYMIALND1wz/xNR9icj+JE9zU9aT
RIar25wa7XZSzF3c9K9wnsJZsJwn8BwnfYw81AxzFXHqx+RpPz0pWgSsXiMxdfa9
04uJqLt5/r1YbF/gJoR+0Cntg9dsWqJTD8YuwGClj3tJ7QSldtI1/isNMziSYRSC
k68Lyy2bvpfvKfDeKdSbHUufs67cXVcLPgxpF1KqgbeN7AWzjxIKlYKxGnDABEVa
KrgtH++piuF1l05cRVXMeLkEpNPxMiwIGJKnu7fpW3IkqAFZ4HvqTm/hyeNf0OnC
0aAve5UnZErVDHVRspkcAQP/8jVBfZU1bscJo361b+XztV2z3p0d3Ut2kTWXMMs=
=aW+/

-----END PGP SIGNATURE-----

## Changelog:

v0.5.0:

	- Fix error with relative import of xml databases

v0.5.0:

	- First release from v0.4.0
	
	- Ported from python2.7 to python3.4
