# burcat.py

# This module extracts the information provided in "Third millenium
# Ideal Gas and Condensed Phase Thermochemical Database for Combustion
# with Updates from Active Thermochemical Tables. A. Burcat and B
# Ruscic. It needs the actual database BURCAT_THR.xml (xml format) to
# run.

# Guillem Borrell i Nogueras
#
# Funding by Vulcano Sadeca

# TODO: Introduce the exception in all try-except clauses

from __future__ import division
try:
    from xml.etree.ElementTree import parse
except ImportError:
    from elementtree import parse
from numpy import empty,array,dot,log
from pkg_resources import Requirement, resource_filename
import doctest
import copy
import os

# Universal gas constant R
R = 8.314472

class Element(object):
    """
    This is a helper class.  It is intended to be created via an
    Elementdb object but you can use it by your own. Take a look at
    Elementdb class.

    Units are K, J, kg... conversion functions are provided in the
    external units module.

    One extra feature not explained in Elementdb documentation is that
    it contains the number of each atom, useful for computing chemical
    reactions.

    >>> db = Elementdb()
    >>> weird = db.getelementdata("C8H6O2")
    >>> print weird.elements
    [('C', 8), ('H', 6), ('O', 2)]
    """
    def __init__(self,formula,Tmin_,_Tmax,mm,hfr,elements):
        self.formula = formula
        self.Tmin_ = Tmin_
        self._Tmax = _Tmax
        self.mm = mm
        self.hfr = hfr
        self.elements = elements

    def density(self,p,T):
        """
        Density.
        """
        return p*self.mm/R/T

    def cpo(self,T):
        """
        Calculates the specific heat capacity in J/mol K
        """
        # I know perfectly that the most efficient way of evaluatin
        # polynomials is recursively but I want the implementation to
        # be as explicit as possible
        Ta = array([1,T,T**2,T**3,T**4],'d')
        if T > 200 and T <= 1000:
            return dot(self.Tmin_[:5],Ta)*R
        elif T >1000 and T < 6000:
            return dot(self._Tmax[:5],Ta)*R
        else:
            raise ValueError("Temperature out of range")

    def cp_(self,T):
        """
        Computes the specific heat capacity in J/kg K for a given temperature
        """
        return self.cpo(T)/self.mm

    @property
    def cp(self):
        """
        Computes the specific heat capacity in J/kg K at 298 K (Reference T)
        """
        return self.cp_(298)

    def ho(self,T):
        """
        Computes the sensible enthalpy in J/mol
        """
        Ta = array([1,T/2,T**2/3,T**3/4,T**4/5,1/T],'d')
        if T > 200 and T < 1000:
            return dot(self.Tmin_[:6],Ta)*R*T
        elif T >1000 and T < 6000:
            return dot(self._Tmax[:6],Ta)*R*T
        else:
            raise ValueError("Temperature out of range")

    def h(self,T):
        """
        Computes the total enthalpy in J/kg
        """
        return self.cp_(T)*T
        
    def so(self,T):
        """
        Computes enthropy in J/mol K
        """
        Ta = array([log(T),T,T**2/2,T**3/3,T**4/4,0,1],'d')
        if T > 200 and T < 1000:
            return dot(self.Tmin_,Ta)*R
        elif T >1000 and T < 6000:
            return dot(self._Tmax,Ta)*R
        else:
            raise ValueError("Temperature out of range")

    def go(self,T):
        """
        Computes the Gibbs free energy from the sensible enthalpy in
        J/mol
        """
        if T > 200 and T < 6000:
            return self.ho(T)-self.so(T)*T
        else:
            raise ValueError("Temperature out of range")
        

    def __repr__(self):
        return """<element> %s"""%(self.formula)

    def __str__(self):
        return """<element> %s"""%(self.formula)

    def __unicode__(self):
        return u"""<element> %s"""%(self.formula)


class Mixture(object):
    """
    Class that models a gas mixture. By now only volume (molar)
    composition is supported.

    You can iterate through all its elements.  The item returned is a
    tuple containing the element and the amount.

    >>> db = Elementdb()
    >>> mix = db.getmixturedata([("O2 REF ELEMENT",20.9476),\
    ("N2  REF ELEMENT",78.084),\
    ("CO2",0.0319),\
    ("AR REF ELEMENT",0.9365),\
    ])
    >>> for e in mix: print e
    (<element> O2 REF ELEMENT, 20.947600000000001)
    (<element> N2  REF ELEMENT, 78.084000000000003)
    (<element> CO2, 0.031899999999999998)
    (<element> AR REF ELEMENT, 0.9365)

    You can get elements either by index or by value.

    >>> print mix['CO2']
    (<element> CO2, 0.031899999999999998)

    You can also delete components of a mixture.  Needed by the
    MoistAir class
    
    >>> mix.delete('CO2')
    >>> print mix
    <Mixture>:
        O2 REF ELEMENT at 20.9476
        N2  REF ELEMENT at 78.084
        AR REF ELEMENT at 0.9365
    """
    def __init__(self,config='vol'):
        self.mix = list()
        self.config = config
        self.idx = 0

    # The following two functions are an iterator. Its purpose is to
    # be able to iterate throug all the elements of a mix.
    def __iter__(self):
        return self

    def next(self):
        try:
            rv = self.mix[self.idx]
            self.idx += 1
            return rv

        except:
            self.idx = 0
            raise StopIteration

    def __getitem__(self,i):
        if type(i) == type(int()):
            return self.mix[i]

        if type(i) == type(str()):
            elem = (None,None)
            for e in self.mix:
                if i == e[0].formula: elem = e

            return elem

    def add(self,component,prop):
        self.mix.append((component,prop))


    def delete(self,formula):
        erased = False
        for e in self.mix:
            if e[0].formula == formula:
                self.mix.remove(e)
                erased = True

        if not erased: raise ValueError("Not a component")

    @property
    def mm(self):
        """
        Computes the equivalent molar mass for a mix
        
        .. math:: 

          M_m = \\frac{1}{N_m} \\sum_i N_i M_i
        """
        if self.config == 'vol':
            Nm = 0
            Mm = 0
            for comp in self.mix:
                Nm += comp[1]

            for comp in self.mix:
                Mm += comp[1]*comp[0].mm

            return Mm/Nm


    def density(self,p,T):
        """
        Computes the density for a given mix of gases

        The equivalent R for a mix is :math:`R_m = \\frac{R_u}{M_n}`,
        where :math:`M_n` is the equivalent molar mass for the mix.
        """
        # TODO: There is a bug in this routine.  Result is not correct.
        return R/self.mm*T/p

    def extensive(self,attr,T):
        """
        Computes the extensive value for a mix.  Remember that an
        extensive value depends on the amount of matter. Enthalpy and
        volume are extensive values.

        .. math::

          ext = \\frac{1}{N_m M_m} \\sum_i N_i M_i ext_i
        """
        if self.config == 'vol':
            Nm = 0
            Mm = 0
            ext = 0
            for comp in self.mix:
                Nm += comp[1]
                
            for comp in self.mix:
                Mm += comp[1]*comp[0].mm

            Mm = Mm/Nm

            for comp in self.mix:
                # Tricky use of getattr function to avoid cutting and
                # pasting several times the very same code
                iattr = getattr(comp[0],attr)
                ext += comp[1] * comp[0].mm * iattr(T)

            return ext/Nm/Mm

    def cp_(self,T):
        """
        Computes the heat capacity at a given temperature

        """
        return self.extensive('cp_',T)

    @property
    def cp(self):
        """
        Computes the heat capacity
        """
        return self.extensive('cp_',298.15)
        
    def ho(self,T):
        return self.extensive('ho',T)

    def h(self,T):
        return self.cp_(T)*T

    def so(self,T):
        return self.extensive('so',T)

    def go(self,T):
        return self.extensive('go',T)

    def __repr__(self):
        str="<Mixture>:"
        for comp in self.mix:
            str +="\n    %s at %s"%(comp[0].formula,comp[1])

        return str

    def __str__(self):
        str="<Mixture>:"
        for comp in self.mix:
            str +="\n    %s at %s"%(comp[0].formula,comp[1])

        return str


    def __unicode__(self):
        str= u"<Mixture>:"
        for comp in self.mix:
            str += u"\n    %s at %s"%(comp[0].formula,comp[1])

        return str


class Elementdb(object):
    """
    Class that reads the Alexander Burcat's thermochemical database
    for combustion.

    >>> db = Elementdb()
    >>> oxygen = db.getelementdata("O2 REF ELEMENT")
    >>> print oxygen
    <element> O2 REF ELEMENT
    >>> print 'molar mass',oxygen.mm
    molar mass 0.0319988
    >>> print 'heat capacity',oxygen.cp
    heat capacity 918.078952423

    The reference temperature for enthalpy is 298.15 K

    >>> print 'enthalpy',oxygen.ho(298.15)
    enthalpy 1.94293914332e-05
    >>> print 'enthropy',oxygen.so(298)
    enthropy 205.133745795
    >>> print 'gibbs free energy',oxygen.go(298)
    gibbs free energy -61134.2629008

    There's a search function.  It is very useful because some names
    are a bit tricky.  Well, not this one.

    >>> db.search("AIR")
    ['AIR']
    >>> air = db.getelementdata("AIR")
    >>> print 'air molar mass',air.mm
    air molar mass 0.02896518
    >>> print 'heat capacity',air.cp
    heat capacity 1004.77625096
    >>> print air.density(101325,298)
    1.1845186553

    The element database can create also mixtures.  It returns an
    instance of Mixture object that can give you the same as the
    Element class for any mixture.

    >>> mix = db.getmixturedata([("O2 REF ELEMENT",20.9476),\
    ("N2  REF ELEMENT",78.084),\
    ("CO2",0.0319),\
    ("AR REF ELEMENT",0.9365),\
    ])
    >>> print mix
    <Mixture>:
        O2 REF ELEMENT at 20.9476
        N2  REF ELEMENT at 78.084
        CO2 at 0.0319
        AR REF ELEMENT at 0.9365
    >>> print mix.cp
    1004.72217065
    >>> print mix.mm
    0.028965116031
    """
    def __init__(self):
        """
        The database file is read when the class is instantiated.
        This is terribly slow as the database is more than 2MB.
        Create the instance and the elements at boot, otherwise be
        prepared to face huge computation times.
        """
        try:
            # try to open the local file, it does not raise an
            # exception on a development environment
            database = open("BURCAT_THR.xml",'r')
        except IOError:
            # Fallback to pkg_resources when thermopy is an installed
            # module
            dbname = resource_filename(
                Requirement.parse("thermopy"),'thermopy/BURCAT_THR.xml')
            database = open(dbname,'r')
            
        tree = parse(database)
        self.db = tree.getroot()

    def search(self,formula):
        """
        List all the species containing a string. Helpful for
        interactive use of the database.
        """
        matches = []
        for specie in self.db:
            try:
                if formula in specie.find("phase").find("formula").text:
                    matches.append(specie.find("phase").find("formula").text)
            except:
                pass

        return matches

    def getelementdata(self,formula):
        """
        Returns an element instance given the name of the element.
        """
        Tmin_ = empty((7),'d')
        _Tmax = empty((7),'d')
        comp = []
        for specie in self.db:
            try:
                if formula == specie.find("phase").find("formula").\
                        text:
                    phase = specie.find("phase")
                    coefficients = phase.find("coefficients")
                    low = coefficients.find("range_Tmin_to_1000") 
                    for (i,c) in zip(range(7),low):
                        Tmin_[i] = float(c.text)

                    high = coefficients.find("range_1000_to_Tmax")
                    for (i,c) in zip(range(7),high):
                        _Tmax[i] = float(c.text)

                    elements = phase.find("elements")
                    elements = elements.getchildren()
                    for elem in elements:
                        it = elem.items()
                        # First is name of element, second is number
                        # of atoms
                        comp.append((it[0][1],int(it[1][1])))

                    mm = float(phase.find("molecular_weight").text)/1000
                    hfr = float(coefficients.find("hf298_div_r").text)
                    
                    return Element(formula,Tmin_,_Tmax,mm,hfr,comp)
                    
            except:
                pass

    def getmixturedata(self,components):
        """
        Creates a mixture of components given a list of tuples
        containing the formula and the volume percent
        """
        mixture = Mixture()
        for comp in components:
            mixture.add(self.getelementdata(comp[0]),comp[1])

        return mixture




if __name__ == '__main__':
    # Move all doctests to py.test
    db = Elementdb()
    mix = db.getmixturedata([("O2 REF ELEMENT",20.9476),
                             ("N2  REF ELEMENT",78.084),
                             ("CO2",0.0319),
                             ("AR REF ELEMENT",0.9365),
                             ("O2 REF ELEMENT",1.2)])
    mix.aggregate()

