from __future__ import division
from burcat import Mixture
from scipy import optimize
from copy import deepcopy

supported = 'C H N O'.split()

# TODO: Check that adiabatic flame temperature is ok

def balance(fuel,am,phi):
    """
    Function that balances the combustion equation given any simple
    fuel element

    :fuel:

      Formula of a Burcat's database fuel as a string. The fuel can
      only contain C, H, O and N.  Other atoms will be ignored.

    :phi:

      Equivalence ratio for air.
    """
    reactants={'fuel':am,'N2':0,'O2':0}
    products={'N2':0,'CO2':0,'H2O':0,'O2':0}
    atoms={'C':0,'H':0,'N':0,'O':0}

    lamb = 1/phi
    for element in fuel.elements:
        if element[0] in supported:
            atoms[element[0]]=element[1]
            
    k = lamb*(atoms['C']+atoms['H']/4+atoms['O']/2)
    reactants['O2'] = am*k
    reactants['N2'] = am*3.76*k

    products['N2'] = am*(atoms['N']+3.76*k)
    products['CO2'] = am*(atoms['C'])
    products['H2O'] = am*(atoms['H']/2)
    products['O2'] = am*((lamb-1)*(k/lamb))

    return (reactants,products)

def balance_mix(fuels,phi):
    """
    function that balances the combustion equation given a mix of
    fuels.

    :fuels:

      type Mixture. Simple fuels formed only by C, H, O and N and the
      amount of each one

    :phi:

      Equivalence ratio for air.
    """
    reactants = {'N2':0,'O2':0}
    products={'N2':0,'CO2':0,'H2O':0,'O2':0}
    atoms={'C':0,'H':0,'N':0,'O':0}

    for fuel in fuels:
        (dreac,dprod)=balance(fuel[0],fuel[1],phi)
        reactants['O2'] += dreac['O2']
        reactants['N2'] += dreac['N2']
        products['N2'] += dprod['N2']
        products['CO2'] += dprod['CO2']
        products['H2O'] += dprod['H2O']
        products['O2'] += dprod['O2']

    return(reactants,products)

class SimpleCombustor(object):
    """
    This class models a simple combustor that uses fuel as a reductor
    and air as a single oxidizer. The combustion is complete, no CO
    nor radicals are formed.

    It only supports fuels with C, N, H and O.  If you put a more
    complicated fuel it will ignore the rest of atoms to balance the
    reaction.

    """
    def __init__(self,fuel,phi,db):
        """
        :fuel:
        
          The formula of the fuel
        
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

        oxygen = db.getelementdata("O2 REF ELEMENT")
        nitrogen = db.getelementdata("N2  REF ELEMENT")
        carbondiox = db.getelementdata("CO2")
        water = db.getelementdata("H2O")
        
        (rdict,pdict) = balance(fuel,1,phi)
        
        self.reactants.add(fuel,rdict['fuel'])
        self.reactants.add(oxygen,rdict['O2'])
        self.reactants.add(nitrogen,rdict['N2'])

        self.products.add(nitrogen,pdict['N2'])
        self.products.add(carbondiox,pdict['CO2'])
        self.products.add(water,pdict['H2O'])
        self.products.add(oxygen,pdict['O2'])

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

    @property
    def lower_heating_value(self):
        return self.heat_of_comb(423.15)

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


class Combustor(object):
    """
    Combustor that is able to characterize the combustion of a mixture
    of fuels
    """
    def __init__(self,fuels,phi,db):
        self.fuels = deepcopy(fuels)
        self.reactants = fuels
        self.products = Mixture()

        oxygen = db.getelementdata("O2 REF ELEMENT")
        nitrogen = db.getelementdata("N2  REF ELEMENT")
        carbondiox = db.getelementdata("CO2")
        water = db.getelementdata("H2O")

        (dreac,dprod) = balance_mix(fuels,phi)

        self.reactants.add(oxygen,dreac['O2'])
        self.reactants.add(nitrogen,dreac['N2'])

        self.products.add(nitrogen,dprod['N2'])
        self.products.add(carbondiox,dprod['CO2'])
        self.products.add(water,dprod['H2O'])
        self.products.add(oxygen,dprod['O2'])

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

        return float(hreac-hprod)/self.fuels.mm

    @property
    def lower_heating_value(self):
        return self.heat_of_comb(423.15)

    def adiabatic_flame_temp(self,T):
        """
        This is the adiabatic flame temp for the given mixtures of
        reactants and products.  If you want the true adiabatic flame
        temperature remember to set the equivalence ratio to 1.
        Otherwise you will always get lower temperatures.
        """
        dh = self.heat_of_comb(T)*self.fuels.mm
        
        f = \
        lambda Tg: self.products[0][1]*self.products[0][0].cpo(Tg)+\
            self.products[1][1]*self.products[1][0].cpo(Tg)+\
            self.products[2][1]*self.products[2][0].cpo(Tg)+\
            self.products[3][1]*self.products[3][0].cpo(Tg)-dh/(Tg-T)
        
        return optimize.fsolve(f,1000)


def test_simplecombustor():
    from burcat import Elementdb
    db = Elementdb()
    methane = db.getelementdata("CH4   RRHO")
    combustor = SimpleCombustor(methane,1.1,db)
    assert combustor.heat_of_comb(298.15) == 50027136.34030433

    # Test Ta
    butane = db.getelementdata('C4H10 n-butane')
    combustor = SimpleCombustor(butane,1,db)
    assert combustor.heat_of_comb(298.15) == 45720359.22491768
    
def test_combustor():
    from burcat import Elementdb
    db = Elementdb()
    fuels = db.getmixturedata([("CH4   RRHO",0.9168),
                               ("C2H6",0.0686),
                               ("C3H8",0.0070),
                               ("C4H10 n-butane",0.0011)])

    combustor = Combustor(fuels,1,db)
    assert combustor.heat_of_comb(423.15) == 49245710.116662093


