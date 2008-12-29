from __future__ import division
from burcat import Mixture
from scipy import optimize
from numpy import float64

class SimpleCombustor(object):
    """
    This class models a simple combustor that uses fuel as a reductor
    and air as a single oxidizer. The combustion is complete.

    It only supports fuels with C, N, H and O.  If you put a more
    complicated fuel it will ignore the rest of atoms to balance the
    reaction.

    """
    supported = 'C H N O'.split()

    def __init__(self,fuel,phi,db):
        """
        :phi:

          Equivalence ratio

        :db:

          Element database.  It is passed as argument because it has
          to be allocated only once.

        >>> from burcat import Elementdb
        >>> db = Elementdb()
        >>> methane = db.getelementdata("CH4   RRHO")
        >>> combustor = SimpleCombustor(methane,1,db)
        >>> print combustor.products.cp
        >>> print combustor.heat_of_comb(298.15)
        50027136.3415
        >>> print combustor.adiabatic_flame_temp(298.15)
        """

        # REACTANTS
        self.reactants = Mixture()
        self.products = Mixture()
        
        self.reactants.append(fuel,1.)
        oxygen = db.getelementdata("O2 REF ELEMENT")
        nitrogen = db.getelementdata("N2  REF ELEMENT")
        
        atoms={'C':0,'H':0,'N':0,'O':0}

        lamb = 1/phi
        for element in fuel.elements:
            if element[0] in self.supported:
                atoms[element[0]]=element[1]
            
        k = lamb*(atoms['C']+atoms['H']/4+atoms['O']/2)
        self.reactants.append(oxygen,k)
        self.reactants.append(nitrogen,3.76*k)

        # PRODUCTS
        self.products = Mixture()
        carbondiox = db.getelementdata("CO2")
        water = db.getelementdata("H2O")
        self.products.append(nitrogen,atoms['N']+3.76*k)
        self.products.append(carbondiox,atoms['C'])
        self.products.append(water,atoms['H']/2)
        self.products.append(oxygen,(lamb-1)*(k/lamb))

    def heat_of_comb(self,T):
        """
        Calculates the heat of combustion per kg of fuel. Checked ok
        """
        hreac = float(0)
        for reac in self.reactants:
            hreac += reac[0].ho(T)*reac[1]
            
        hprod = float(0)
        for prod in self.products:
            hprod += prod[0].ho(T)*prod[1]

        return float(hreac-hprod)/self.reactants[0][0].mm


    def adiabatic_flame_temp(self,T):
        """
        This is the adiabatic flame temp for the given mixtures of
        reactants and products.  If you want the true adiabatic flame
        temperature remember to set the equivalence ratio to 1.
        Otherwise you will always get lower temperatures.
        """
        dh = self.heat_of_comb(T)*self.reactants[0][0].mm
        
        f = \
        lambda Tg: self.products[0][1]*self.products[0][0].cpo(Tg)+\
            self.products[1][1]*self.products[1][0].cpo(Tg)+\
            self.products[2][1]*self.products[2][0].cpo(Tg)+\
            self.products[3][1]*self.products[3][0].cpo(Tg)-dh/(Tg-T)
        
        return optimize.fsolve(f,1000)


class CustomCombustor(object):
    """
    Calculates the heat of combustion per kg of fuel. Checked ok
    """
    def __init__(self,reactants,products):
        pass

def run():
    from burcat import Elementdb
    db = Elementdb()
    methane = db.getelementdata("CH4   RRHO")
    combustor = SimpleCombustor(methane,1.1,db)
    print 'heat of combustion',combustor.heat_of_comb(298.15)

    # Test Ta
    butane = db.getelementdata('C4H10 n-butane')
    combustor = SimpleCombustor(butane,1,db)
    print 'Ta',combustor.adiabatic_flame_temp(298.15)
    



            
if __name__ == "__main__":
    #from doctest import testmod
    #testmod()
    run()

    
