u"""
Units conversion module.

Classes:

    Temperature: temperature object.

    Pressure: pressure object.

    Enthalpy: enthalpy object.

    Length: length object.

    Massflow: mass flow object.

    Massflowrate: massflowrate object.

    Energy: energy object.

"""
import thermopy.constants as constants


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
    37.777777777777...
    >>> T.F
    99.999999999999...

    You can compute with temperatures because inherits from the float
    built-in

    >>> T1 = Temperature(200)
    >>> T2 = Temperature(0).unit('C')
    >>> round(T1+T2, 2)
    473.15

    If you don't want to use the class' attribute you can use the
    function `getattr` to get a value using the unit code.

    >>> getattr(T,'C')
    37.77777777777...
    """

    def __init__(self, data):
        u"""Initialize temperature object."""
        self.data = float(data)

    @classmethod
    def __factory(cls, data):
        """
        This factory makes that any returned value is a measure
        instead of a float.
        """
        return cls(data)

    def unit(self, units='K'):
        u"""Set unit for temperature."""
        if units == 'K':
            return self.__factory(self.data)
        elif units == 'C':
            return self.__factory(constants.C2K(self.data))
        elif units == 'F':
            return self.__factory(constants.F2K(self.data))
        else:
            raise ValueError("Wrong temperature input code")

    @property
    def C(self):
        u"""Property of Celsius temperature unit."""
        return self.__factory(constants.K2C(self.data))

    @property
    def F(self):
        u"""Property of Fahrenheit temperature unit."""
        return self.__factory(constants.K2F(self.data))


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

    def __init__(self, data):
        u"""Initialize pressure object."""
        self.data = float(data)

    @classmethod
    def __factory(cls, data):
        """
        This factory makes that any returned value is a Measure
        instead of a float.
        """
        return cls(data)

    def unit(self, units='Pa'):
        u"""Set unit for pressure."""
        if units == 'Pa':
            return self.__factory(self.data)
        elif units == 'MPa':
            return self.__factory(constants.mega * self.data)
        elif units == 'bar':
            return self.__factory(self.data * constants.bar)
        elif units == 'psi':
            return self.__factory(self.data * constants.psi)
        elif units == 'atm':
            return self.__factory(self.data * constants.atm)
        elif units == 'mmwc':
            return self.__factory(self.data * (constants.torr * 1000 / 13534))
        elif units == 'torr':
            return self.__factory(self.data * constants.torr)
        else:
            raise ValueError("wrong pressure unit input code")

    @property
    def MPa(self):
        u"""Property of MPa pressure unit."""
        return self.__factory(self.data / constants.mega)

    @property
    def bar(self):
        u"""Property of bar pressure unit."""
        return self.__factory(self.data / constants.bar)

    @property
    def psi(self):
        u"""Property of psi pressure unit."""
        return self.__factory(self.data / constants.psi)

    @property
    def atm(self):
        u"""Property of atm pressure unit."""
        return self.__factory(self.data / constants.atm)

    @property
    def mmwc(self):
        u"""Property of mmwc pressure unit."""
        return self.__factory(self.data / (constants.torr * 1000 / 13534))

    @property
    def torr(self):
        u"""Property of torr pressure unit."""
        return self.__factory(self.data / constants.torr)

# HUNITS = ['si', 'kJkg', 'kcalkg', 'Btulb']


class Enthalpy(float):
    """
    Class that models an enthalpy measure with conversion utilities

    Supported units are

    * Joule per kg (default)

    * Kilojoule per kg (kJkg)

    * Kilocalorie per kg (kcalkg)

    * BTU per pound (Btulb)

    >>> h = Enthalpy(1000)
    >>> h.kJkg
    1.0
    >>> h.kcalkg
    0.2390057361376...
    >>> h.Btulb
    0.42992261392949266
    """

    def __init__(self, data):
        u"""Initialize enthalpy object."""
        self.data = float(data)

    @classmethod
    def __factory(cls, data):
        """
        This factory makes that any returned value is a measure
        instead of a float.
        """
        return cls(data)

    def unit(self, units='si'):
        u"""Set unit for Enthalpy."""
        if units == 'si':
            return self.__factory(self.data)
        elif units == 'kJkg':
            return self.__factory(self.data * constants.kilo)
        elif units == 'kcalkg':
            return self.__factory(self.data * constants.calorie *
                                  constants.kilo)
        elif units == 'Btulb':
            return self.__factory(self.data * constants.Btu / constants.lb)
        raise ValueError("wrong enthalpy unit input code")

    @property
    def kJkg(self):
        u"""Property of KJkg pressure unit."""
        return self.__factory(self.data / constants.kilo)

    @property
    def kcalkg(self):
        u"""Property of kcalkg pressure unit."""
        return self.__factory(self.data / constants.kilo / constants.calorie)

    @property
    def Btulb(self):
        u"""Property of Btulb pressure unit."""
        return self.__factory(self.data * constants.lb / constants.Btu)


class Length(float):

    """
    Class that models a length measure with conversion utilities

    Supported units are

    * meter (default)

    * millimeter (mm)

    * inch (inch)

    * foot (ft)

    >>> l = Length(1).unit('inch')
    >>> round(l.mm, 1)
    25.4
    >>> l.ft
    0.0833333333333...
    >>> round(l, 4)
    0.0254
    """

    def __init__(self, data):
        u"""Initialize length object."""
        self.data = float(data)

    @classmethod
    def __factory(cls, data):
        """
        This factory makes that any returned value is a measure
        instead of a float.
        """
        return cls(data)

    def unit(self, units='m'):
        u"""Set unit for length."""
        if units == 'm':
            return self.__factory(self.data)
        elif units == 'mm':
            return self.__factory(self.data * constants.milli)
        elif units == 'inch':
            return self.__factory(self.data * constants.inch)
        elif units == 'ft':
            return self.__factory(self.data * constants.foot)
        else:
            raise ValueError("wrong length unit input code")

    @property
    def mm(self):
        u"""Property of mm unit."""
        return self.__factory(self.data / constants.milli)

    @property
    def inch(self):
        u"""Property of inch unit."""
        return self.__factory(self.data / constants.inch)

    @property
    def ft(self):
        u"""Property of feet unit."""
        return self.__factory(self.data / constants.foot)


class Massflow(float):
    """
    Class that models a mass flow measure with conversion utilities

    Supported units are

    * kg per second (default)

    * kg per hour (kgh)

    * pounds per second (lbs)

    * pounds per hour (lbh)
    """

    def __init__(self, data):
        u"""Initialize massflow object."""
        self.data = float(data)

    @classmethod
    def __factory(cls, data):
        """
        This factory makes that any returned value is a measure
        instead of a float.
        """
        return cls(data)

    def unit(self, units='kgs'):
        u"""Set unit for massflow."""
        if units == 'kgs':
            return self.__factory(self.data)
        elif units == 'kgh':
            return self.__factory(self.data / constants.hour)
        elif units == 'lbs':
            return self.__factory(self.data * constants.lb)
        elif units == 'lbh':
            return self.__factory(self.data * constants.lb / constants.hour)
        else:
            raise ValueError("wrong massflow unit input code")

    @property
    def kgh(self):
        u"""Property of kgh unit."""
        return self.__factory(self.data * constants.hour)

    @property
    def lbs(self):
        u"""Property of lbs unit."""
        return self.__factory(self.data / constants.lb)

    @property
    def lbh(self):
        u"""Property of lbh unit."""
        return self.__factory(self.data * constants.hour / constants.lb)


class Massflowrate(float):
    """
    Class that models a mass flow measure with conversion utilities

    Supported units are

    * :math:`\\frac{kg}{s\ m^2}` (default)

    * :math:`\\frac{lb}{s\ ft^2}` (Btu)
    """

    def __init__(self, data):
        u"""Initialize massflowrate object."""
        self.data = float(data)

    @classmethod
    def __factory(cls, data):
        """
        This factory makes that any returned value is a measure
        instead of a float.
        """
        return cls(data)

    def unit(self, units='default'):
        u"""Set unit for massflowrate."""
        if units == 'default':
            return self.__factory(self.data)
        elif units == 'Btu':
            return self.__factory(self.data * constants.lb /
                                  constants.foot ** 2)
        else:
            raise ValueError("wrong massflow unit input code")

    @property
    def Btu(self):
        u"""Property of Btu unit."""
        return self.__factory(self.data * constants.foot ** 2 / constants.lb)


class Energy(float):
    """Energy in J, Btu, cal, kWh"""
    def __init__(self, data):
        u"""Initialize energy object."""
        self.data = float(data)

    @classmethod
    def __factory(cls, data):
        """
        This factory makes that any returned value is a measure
        instead of a float.
        """
        return cls(data)

    def unit(self, units='J'):
        u"""Set unit for energy."""
        if units.upper() == 'BTU':
            return self.__factory(self.data * constants.Btu)
        elif units == 'cal':
            return self.__factory(self.data * constants.calorie)
        elif units == 'kWh':
            return self.__factory(self.data * constants.kWh)

    @property
    def Btu(self):
        u"""Property of Btu unit."""
        return self.__factory(self.data / constants.Btu)

    @property
    def cal(self):
        u"""Property of cal unit."""
        return self.__factory(self.data / constants.calorie)

    @property
    def kWh(self):
        u"""Property of KWh unit."""
        return self.__factory(self.data / constants.kWh)
