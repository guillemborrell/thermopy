from __future__ import division
from iapws import Water
from units import Tu, pu, hu


# Water at triple point
pt = 612
Tt = 273.16

class MoistAir(object):
    """
    Class that models a moist gas.  The computations in this class are
    a bit tricky because the enthalpy computations for steam far from
    atmospheric pressure using the bcura data have severe
    deviations. Considering water an ideal gas is a too strong
    assumption. This means that water data is computed using the IAPWS
    data for water an steam.

    The trickiest part comes with the fact that the enthalpy reference
    for the IAPWS tables is the water's triple point when the enthalpy
    reference for the bcura tables is the absolute zero.

    The IAPWS reference is the leading one.

    MoistAir does not inherit from Mixture because of this.
    """
    def __init__(self,gas):
        self.moist = Water()

        (self.water,self.qwater) = gas['H2O']
        if self.water == None: raise ValueError("No water in this gas")

        gas.delete('H2O')
        self.gas = gas

        # amount of gas (in case it is not 100-qwater)
        self.qgas = 0
        for e in self.gas:
            self.qgas += e[1]

        self.w = self.water.mm/self.gas.mm*self.qwater/self.qgas
        
    def phi(self,p,T):
        """
        Relative moisture given pressure and temperature.
        """
        ya = self.qgas/(self.qgas+self.qwater)
        return self.gas.mm*ya*p*self.w/(self.water.mm*self.moist.psat(T))


    def wet_bulb_T(self,p):
        """
        Wet bulb temperature for a given pressure
        """
        yw = self.qwater/(self.qgas+self.qwater)
        return self.moist.Tsat(yw*p)

    def h(self,p,T):
        raise NotImplementedError

    def __repr__(self):
        return "<Moist Gas>:\n  Gas:\n"+self.gas.__repr__()


    def __unicode__(self):
        return u"<Moist Gas>:\n  Gas:\n"+self.gas.__unicode__()


def test_wark():
    """
    This function runs the 10.7 example from Wark and Richard's
    Thermodynamics, spanish translation.
    """
    from burcat import Elementdb
    from units import pu,Tu
    db = Elementdb()
    gas = db.getmixturedata([("AIR",1),("H2O",0.015)])
    ma = MoistAir(gas)
    assert ma.w == 0.0093294500500255822
    assert ma.phi(pu(14.7,f='psi'),Tu(70,f='F')) == 0.59790008408358786
    assert ma.wet_bulb_T(pu(14.7,f='psi')) == 286.14757997335232
    

