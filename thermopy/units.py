from constants import *

class Temperature(float):
    """
    Class that models a temperature measure with conversion utilities

    Supported units are

    * Kelvin

    * Celsius

    * Fahrenheit

    Normal instantiation is a temperature in Kelvin

    >>> T = Temperature(100)
    >>> T
    100.0

    But you can instantiate and specify if unit is Celsius or
    Fahrenheit

    >>> T = Temperature(100).unit('F')
    >>> T
    310.92777777777775

    Unit conversion is as easy as it gets.

    >>> T.C
    37.777777777777771
    >>> T.F
    99.999999999999986

    You can compute with temperatures because inherits from the float
    built-in

    >>> T1 = Temperature(200)
    >>> T2 = Temperature(0).unit('C')
    >>> T1+T2
    473.14999999999998
    """
    def __init__(self,data):
        float.__init__(self,float(data))
        self.data = float(data)

    @classmethod
    def factory(cls,data):
        """
        This factory makes that any returned value is a Temperature
        instead of a float.
        """
        return cls(data)

    def unit(self,units='K'):
        if units == 'K':
            return self.factory(self.data)
        elif units == 'C':
            return self.factory(C2K(self.data))
        elif units == 'F':
            return self.factory(F2K(self.data))
        else:
            raise ValueError("Wrong temperature input code")
        
    @property
    def C(self):
        return self.factory(K2C(self.data))

    @property
    def F(self):
        return self.factory(K2F(self.data))


class Pressure(float):
    """
    Class that models a Pressure measure with conversion utilities

    Suported units are

    * Pascal (Pa)

    * Mega Pascal (MPa)

    * Bar (bar)

    * Pound per square inch (psi)

    * Atmosphere (atm)

    * Millimeters of water column (mmwc)

    * Torricelli (torr)

    Normal instantiation is pressure in Pa. How much is an athmosphere?

    >>> p = Pressure(1.0).unit('atm')
    >>> p
    101325.0
    >>> p.torr
    760.0
    >>> p.mmwc
    10285.839999999998
    >>> p.psi
    14.69594877551345
    """
    def __init__(self,data):
        float.__init__(self,float(data))
        self.data = float(data)

    @classmethod
    def factory(cls,data):
        """
        This factory makes that any returned value is a Temperature
        instead of a float.
        """
        return cls(data)

    def unit(self,units='Pa'):
        if units == 'Pa':
            return self.factory(self.data)
        elif units == 'MPa':
            return self.factory(mega*self.data)
        elif units == 'bar':
            return self.factory(self.data*bar)
        elif units == 'psi':
            return self.factory(self.data*psi)
        elif units == 'atm':
            return self.factory(self.data*atm)
        elif units == 'mmwc':
            return self.factory(self.data*(torr*1000/13534))
        elif units == 'torr':
            return self.factory(self.data*torr)
        else:
            raise ValueError("wrong pressure unit input code")

    @property
    def MPa(self):
        return self.factory(self.data/mega)

    @property
    def bar(self):
        return self.factory(self.data/bar)

    @property
    def psi(self):
        return self.factory(self.data/psi)

    @property
    def bar(self):
        return self.factory(self.data/bar)

    @property
    def atm(self):
        return self.factory(self.data/atm)

    @property
    def mmwc(self):
        return self.factory(self.data/(torr*1000/13534))

    @property
    def torr(self):
        return self.factory(self.data/torr)

HUNITS=['si','kJkg','kcalkg','Btulb']

def hu(h,f='si',t='si'):
    """
    Helper function to change enthalpy given units

    >>> hu(970.7,'Btulb','kcalkg')
    539.63867112810715
    """
    if f == 'si':
        h = h

    elif f == 'kJkg':
        h = h*kilo

    elif f == 'kcalkg':
        h = h*calorie*kilo

    elif f == 'Btulb':
        h = h*Btu/lb

    else:
        raise ValueError("wrong enthalpy unit input code")

    if t == 'si':
        return h

    elif t == 'kJkg':
        return h/kilo

    elif t == 'kcalkg':
        return h/calorie/kilo

    elif f == 'Btulb':
        return h*lb/Btu

    else:
        raise ValueError("wrong enthalpy unit output code")

LUNITS=['m','mm','in','ft']

def lu(l,f='m',t='m'):
    """
    Utility function to change length units

    >>> lu(8,f='ft',t='mm')
    2438.3999999999996
    >>> lu(1.5,f='in',t='mm')
    38.099999999999994
    """
    if f == 'm':
        l = l

    elif f == 'mm':
        l = milli*l

    elif f == 'in':
        l = l*inch

    elif f == 'ft':
        l = l*foot

    else:
        raise ValueError("wrong length unit input code")

    if t == 'm':
        return l

    elif t == 'mm':
        return l/milli

    elif t == 'in':
        return l/inch

    elif t == 'ft':
        return l/foot

    else:
        raise ValueError("wrong length unit output code")

MFUNITS=['si','kgh','lbs','lbh']

def mfu(massflow,f='si',t='si'):
    """
    Utility function to change the massflow units.
    
    >>> mfu(1.47,'lbs','si')
    0.66678078389999995
    """
    if f == 'si':
        massflow = massflow

    elif f == 'kgh':
        massflow = massflow/hour

    elif f == 'lbs':
        massflow = massflow*lb

    elif f == 'lbh':
        massflow = massflow*lb/hour

    else:
        raise ValueError("wrong massflow unit input code")

    if t == 'si':
        return massflow

    elif t == 'kgh':
        return massflow*hour

    elif t == 'lbs':
        return massflow/lb

    elif t == 'lbh':
        return massflow*hour/lb

    else:
        raise ValueError("wrong massflow unit output code")

MFRUNITS=['si','btu']

def mfru(massflowrate,f='si',t='si'):
    """
    Utility function to change the massflowrate units.
    """
    if f == 'si':
        massflowrate = massflowrate

    elif f == 'btu':
        massflowrate = massflowrate*lb/foot**2

    else:
        raise ValueError("wrong massflow rate unit input code")

    if t == 'si':
        return massflowrate

    elif t == 'btu':
        return massflowrate*foot**2/lb

    else:
        raise ValueError("wrong massflow rate output code")

DENSITYU=['si']

def densityu(density,f='si',t='si'):
    if f == 'si':
        density = density
    else:
        raise NotImplementedError("No other units supported yet")

    if t == 'si':
        return density
    else:
        raise NotImplementedError("No other units supported yet")


def test_doctest():
    import doctest
    doctest.testmod()


