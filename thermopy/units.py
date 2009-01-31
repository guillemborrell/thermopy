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

    If you don't want to use the class' attribute you can use the
    function `getattr` to get a value using the unit code.

    >>> getattr(T,'C')
    37.777777777777771
    """
    def __init__(self,data):
        float.__init__(self,float(data))
        self.data = float(data)

    @classmethod
    def __factory(cls,data):
        """
        This factory makes that any returned value is a measure
        instead of a float.
        """
        return cls(data)

    def unit(self,units='K'):
        if units == 'K':
            return self.__factory(self.data)
        elif units == 'C':
            return self.__factory(C2K(self.data))
        elif units == 'F':
            return self.__factory(F2K(self.data))
        else:
            raise ValueError("Wrong temperature input code")

    @property
    def K(self):
        return self.__factory(self.data)
        
    @property
    def C(self):
        return self.__factory(K2C(self.data))

    @property
    def F(self):
        return self.__factory(K2F(self.data))


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
    def __factory(cls,data):
        """
        This factory makes that any returned value is a Measure
        instead of a float.
        """
        return cls(data)

    def unit(self,units='Pa'):
        if units == 'Pa':
            return self.__factory(self.data)
        elif units == 'MPa':
            return self.__factory(mega*self.data)
        elif units == 'bar':
            return self.__factory(self.data*bar)
        elif units == 'psi':
            return self.__factory(self.data*psi)
        elif units == 'atm':
            return self.__factory(self.data*atm)
        elif units == 'mmwc':
            return self.__factory(self.data*(torr*1000/13534))
        elif units == 'torr':
            return self.__factory(self.data*torr)
        else:
            raise ValueError("wrong pressure unit input code")


    @property
    def Pa(self):
        return self.__factory(self.data)

    @property
    def MPa(self):
        return self.__factory(self.data/mega)

    @property
    def bar(self):
        return self.__factory(self.data/bar)

    @property
    def psi(self):
        return self.__factory(self.data/psi)

    @property
    def atm(self):
        return self.__factory(self.data/atm)

    @property
    def mmwc(self):
        return self.__factory(self.data/(torr*1000/13534))

    @property
    def torr(self):
        return self.__factory(self.data/torr)

HUNITS=['Jkg','kJkg','kcalkg','Btulb']

class Enthalpy(float):
    """
    Class that models an enthalpy measure with conversion utilities

    Supported units are

    * Joule per kg (Jkg)

    * Kilojoule per kg (kJkg)

    * Kilocalorie per kg (kcalkg)

    * BTU per pound (Btulb)

    >>> h = Enthalpy(1000)
    >>> h.kJkg
    1.0
    >>> h.kcalkg
    0.23900573613766729
    >>> h.Btulb
    0.42992261392949266
    """

    def __init__(self,data):
        float.__init__(self,float(data))
        self.data = float(data)

    @classmethod
    def __factory(cls,data):
        """
        This factory makes that any returned value is a measure
        instead of a float.
        """
        return cls(data)

    def unit(self,units='Jkg'):
        if units == 'Jkg':
            return self.__factory(self.data)
        elif units == 'kJkg':
            return self.__factory(self.data*kilo)
        elif units == 'kcalkg':
            return self.__factory(self.data*calorie*kilo)
        elif units == 'Btulb':
            return self.__factory(self.data*Btu/lb)
        raise ValueError("wrong enthalpy unit input code")

    @property
    def Jkg(self):
        return self.__factory(self.data)

    @property
    def kJkg(self):
        return self.__factory(self.data/kilo)

    @property
    def kcalkg(self):
        return self.__factory(self.data/kilo/calorie)

    @property
    def Btulb(self):
        return self.__factory(self.data*lb/Btu)


class Length(float):
    """
    Class that models a length measure with conversion utilities

    Supported units are

    * meter (m)

    * millimeter (mm)

    * inch (inch)

    * foot (ft)

    >>> l = Length(1).unit('inch')
    >>> l.mm
    25.399999999999999
    >>> l.ft
    0.083333333333333343
    >>> l
    0.025399999999999999
    """
    def __init__(self,data):
        float.__init__(self,float(data))
        self.data = float(data)

    @classmethod
    def __factory(cls,data):
        """
        This factory makes that any returned value is a measure
        instead of a float.
        """
        return cls(data)

    def unit(self,units='m'):
        if units == 'm':
            return self.__factory(self.data)
        elif units == 'mm':
            return self.__factory(self.data*milli)
        elif units == 'inch':
            return self.__factory(self.data*inch)
        elif units == 'ft':
            return self.__factory(self.data*foot)
        else:
            raise ValueError("wrong length unit input code")

    @property
    def m(self):
        return self.__factory(self.data)

    @property
    def mm(self):
        return self.__factory(self.data/milli)

    @property
    def inch(self):
        return self.__factory(self.data/inch)

    @property
    def ft(self):
        return self.__factory(self.data/foot)
    

class Massflow(float):
    """
    Class that models a mass flow measure with conversion utilities

    Supported units are

    * kg per second (kgs)

    * kg per hour (kgh)

    * pounds per second (lbs)

    * pounds per hour (lbh)
    """
    def __init__(self,data):
        float.__init__(self,float(data))
        self.data = float(data)

    @classmethod
    def __factory(cls,data):
        """
        This factory makes that any returned value is a measure
        instead of a float.
        """
        return cls(data)

    def unit(self,units='kgs'):
        if units == 'kgs':
            return self.__factory(self.data)
        elif units == 'kgh':
            return self.__factory(self.data/hour)
        elif units == 'lbs':
            return self.__factory(self.data*lb)
        elif units == 'lbh':
            return self.__factory(self.data*lb/hour)
        else:
            raise ValueError("wrong massflow unit input code")

    @property
    def kgs(self):
        return self.__factory(self.data)

    @property
    def kgh(self):
        return self.__factory(self.data*hour)

    @property
    def lbs(self):
        return self.__factory(self.data/lb)

    @property
    def lbh(self):
        return self.__factory(self.data*hour/lb)
    

class Massflowrate(float):
    """
    Class that models a mass flow measure with conversion utilities

    Supported units are

    * :math:`\\frac{kg}{s\ m^2}` (si)

    * :math:`\\frac{lb}{s\ ft^2}` (Btu)
    """
    def __init__(self,data):
        float.__init__(self,float(data))
        self.data = float(data)

    @classmethod
    def __factory(cls,data):
        """
        This factory makes that any returned value is a measure
        instead of a float.
        """
        return cls(data)

    def unit(self,units='si'):
        if units == 'si':
            return self.__factory(self.data)
        elif units == 'Btu':
            return self.__factory(self.data*lb/foot**2)
        else:
            raise ValueError("wrong massflow unit input code")

    @property
    def si(self):
        return self.__factory(self.data)

    @property
    def Btu(self):
        return self.__factory(self.data*foot**2/lb)


def test_doctest():
    import doctest
    doctest.testmod()


