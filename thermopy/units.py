from constants import *


class Temperature(float):
    def __init__(self,data):
        if unit == 'K':
            float.__init__(self,data)
        elif unit == 'C':
            float.__init__(self,C2K(data))
        elif unit == 'F':
            float.__init__(self,F2K(data))
        else:
            raise ValueError("wrong temperature unit input code")

    @property
    def C(self):
        return K2C(self.data)

    @property
    def F(self):
        return K2F(self.data)

def test_temperature():
    t = Temperature(100)
    assert type(t) == 1
    assert t.F == 1


def Tu(T,f='K',t='K'):
    """
    Helper function to change temperature given units
    """
    if f == 'K':
        T = T
        
    elif f == 'C':
        T = C2K(T)

    elif f == 'F':
        T = F2K(T)

    else:
        raise ValueError("wrong temperature unit input code")

    if t == 'K':
        return T

    elif t == 'C':
        return K2C(T)

    elif t == 'F':
        return K2F(T)

    else:
        raise ValueError("wrong temperature unit output code")

PUNITS=['Pa','MPa','bar','psi','atm','mmwc']

def pu(p,f='Pa',t='Pa'):
    """
    Helper function to change pressure given units
    """
    if f == 'Pa':
        p = p

    elif f == 'MPa':
        p = mega*p

    elif f == 'bar':
        p = p*bar

    elif f == 'psi':
        p = p*psi
        
    elif f == 'atm':
        p = p*atm

    elif f == 'mmwc':
        p = p*(torr*1000/13534)
    
    else:
        raise ValueError("wrong pressure unit input code")

    if t == 'Pa':
        return p

    elif t == 'MPa':
        return p/mega

    elif t == 'bar':
        return p/bar

    elif t == 'psi':
        return p/psi

    elif t == 'atm':
        return p/atm

    elif t == 'mmwc':
        return p/(torr*1000/13534)

    else:
        raise ValueError("wrong pressure unit output code")

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


def test():
    import doctest
    doctest.testmod()


if __name__ == '__main__':
    test()
