"""
Burcat module is intended to give access to the Burcat's Database [1].

This database contains coefficients for thousands of chemicals to be used in
polynomial form. Therefore it allows the user to calculate different
thermodynamic properties for various compounds. Some of these properties are: 1.
Specific heat capacity at constant pressure 2. Enthalpy 3. Entropy 4. Gibbs
energy all of them as functions of temperature for a given compound.

Intended audience from [1]: The database is used by scientists, educators,
engineers and students at all levels, dealing primarily with combustion and air
pollution, jet engines, rocket propulsion, fireworks, but also by researchers
involved in upper atmosphere kinetics, astrophysics, abrasion metallurgy, etc.

This file is deprecated and shall not be maintained. Use the nasa9polynomials
instead.
"""

import os
from xml.etree.ElementTree import parse
import numpy as np
import thermopy.units as units
from thermopy.constants import ideal_gas_constant
_R = ideal_gas_constant[0]


class Compound(object):
    u"""
    Create chemical compounds.

    It is intended to be created via an Elementdb object but you can use it by
    your own. Take a look at Elementdb class.

    Attributes:
        :param cas: (?type?): CAS number of the compound.
        description (??): ?
        formula (??): ?
        aggr_state (??): ?
        mm (??): molar mass of the compound.
        h_formation (): heat of formation.

    Units are in SI on a molar basis.
    """

    def __init__(self, cas, description, reference, formula, elements,
                 aggr_state, T_limit_low, T_limit_high, calc_quality,
                 mm, low_coefs, high_coefs, h_formation):
        """
        Create a chemical compound usually from the set_compound method.

        Args:
            cas (type):
            description ():
            reference ():


        """
        self.cas = cas
        self.description = description
        self._reference = reference
        self.formula = formula
        self._elements = elements
        self.aggr_state = aggr_state
        self._T_limit_low = T_limit_low
        self._T_limit_high = T_limit_high
        self._calc_quality = calc_quality
        self.mm = mm
        self._low_coefs = low_coefs
        self._high_coefs = high_coefs
        self.h_formation = h_formation

    def density_ideal_gas(self, p, T):
        """
        Density in kg/m³.
        """
        return units.Pressure(p) * self.mm / _R / units.Temperature(T)

    def heat_capacity(self, T):
        """
        Calculates the specific heat capacity in J/(mol K).
        """
        Ta = np.array([1, T, T ** 2, T ** 3, T ** 4], dtype='d')
        if T >= self._T_limit_low and T <= 1000:
            return np.dot(self._low_coefs[:5], Ta) * _R
        elif T > 1000 and T <= self._T_limit_high:
            return np.dot(self._high_coefs[:5], Ta) * _R
        else:
            raise ValueError("Temperature out of range")

    def heat_capacity_massic(self, T):
        """
        Computes the specific heat capacity in J/(kg K)
        for a given temperature.
        """
        return self.cpo(T) / self.mm

    def enthalpy(self, T):
        """
        Computes the sensible enthalpy in J/mol.
        """
        Ta = np.array([1, T / 2, T ** 2 / 3, T ** 3 / 4, T ** 4 / 5, 1 / T],
                      'd')
        if T >= self._T_limit_low and T <= 1000:
            partial = np.dot(self._low_coefs[:6], Ta) * _R * T
        elif T > 1000 and T <= self._T_limit_high:
            partial = np.dot(self._high_coefs[:6], Ta) * _R * T
        else:
            raise ValueError("Temperature out of range")

        return partial - self.h_formation

    def enthalpy_massic(self, T):
        """
        Computes the sensible enthalpy in J/kg.
        """
        Ta = np.array([1, T / 2, T ** 2 / 3, T ** 3 / 4, T ** 4 / 5, 1 / T],
                      'd')
        if T >= self._T_limit_low and T <= 1000:
            partial = np.dot(self._low_coefs[:6], Ta) * _R * T / self.mm
        elif T > 1000 and T <= self._T_limit_high:
            partial = np.dot(self._high_coefs[:6], Ta) * _R * T / self.mm
        else:
            raise ValueError("Temperature out of range")

        return partial - self.h_formation

    def enthalpy_engineering(self, T):
        """
        Computes the total enthalpy in J/mol: h = h_formation +
        integral cp(T) dT.\nUseful for heating reacting systems such as:
        MgO + Cl2 -> MgCl2 + 0.5 O2
        from Tref to 500 K.\ndeltaH total = (sum (h_form product) -
        sum (h_form react)) + sum ( enthalpy(500) * prod)_each_product
        """
        return self.h_formation + self.enthalpy(T)

    def entropy(self, T):
        """
        Computes enthropy in J/mol K.
        """
        Ta = np.array([np.log(T), T, T ** 2 / 2, T ** 3 / 3, T ** 4 / 4, 0, 1],
                      'd')
        # right
        if T >= self._T_limit_low and T <= 1000:
            return np.dot(self._low_coefs, Ta) * _R
        elif T > 1000 and T <= self._T_limit_high:
            return np.dot(self._high_coefs, Ta) * _R
        else:
            raise ValueError("Temperature out of range")

    def gibbs_energy(self, T):
        """
        Computes the Gibbs free energy from the sensible enthalpy in
        J/mol.
        """
        if T >= self._T_limit_low and T < self._T_limit_high:
            return self.enthalpy(T) - self.entropy(T) * T
        else:
            raise ValueError("Temperature out of range")

    def __repr__(self):
        return """<element> %s""" % (self.formula)

    def __str__(self):
        return """<element> %s""" % (self.formula)

    def __unicode__(self):
        return u"""<element> %s""" % (self.formula)


class Database(object):
    """
    Class that reads the Alexander Burcat's thermochemical database
    for combustion.
    """

    def __init__(self):
        """
        The database file is read when the class is instantiated.
        This is terribly slow as the database is more than 2MB.
        Create the instance and the elements at boot, otherwise be
        prepared to face huge computation times.
        """
        self.db = parse(str(os.path.dirname(os.path.dirname(__file__)) +
                        '/databases/burcat_thr.xml')).getroot()

    def list_compound(self, cas_or_formula):
        """
        Takes a string or a CAS number as input and output all
        matching results. Helpful for interactive use of the database.
        """
        # determines if it is cas or not
        for char in cas_or_formula:
            if char.isalpha() is True:
                formula = cas_or_formula
                cas = None
                break
        else:
            cas = cas_or_formula
            formula = None
        matches = []
        if cas is not None:
            for specie in self.db.findall('specie'):
                if cas == str(specie.get('CAS')):
                    try:
                        specie.find('phase')
                        for phases in specie.findall('phase'):
                            matches.append(phases.find('formula').text)
                    except:
                        pass
            return matches
        elif formula is not None:
            formula = formula.upper()
            for specie in self.db.findall('specie'):
                try:
                    specie.find('phase')
                    for phases in specie.findall('phase'):
                        if formula in phases.find('formula').text.upper():
                            matches.append(phases.find('formula').text)
                except:
                    pass
            return matches

    def set_compound(self, formula):
        """
        Returns an Element instance given the name of the element.
        """
        formula = formula.upper()
        for specie in self.db.findall('specie'):
            try:
                specie.findall('phase')
                for each_phase in specie.findall('phase'):
                    try:
                        each_phase.findall('formula')
                        for each_formula in each_phase.findall('formula'):
                            if formula == each_formula.text.upper():
                                cas = str(specie.get('CAS'))
                                try:
                                    specie.find('formula_name_structure').find(
                                        'formula_name_structure_1')
                                    description = str(
                                        specie.find(
                                            'formula_name_structure').find(
                                                'formula_name_structure_1'
                                            ).text)
                                except:
                                    description = None
                                try:
                                    specie.find('formula_name_structure_1')
                                    reference = str(specie.find(
                                        'formula_name_structure_1').text)
                                except:
                                    reference = None
                                    elements = []
                                    for elem in each_phase.find('elements'):
                                        elements.append((elem.get('name'),
                                                         int(elem.get(
                                                             'num_of_atoms'))))
                                    aggr_state = str(each_phase.find(
                                        'phase').text)
                                    T_limit_low = float(each_phase.find(
                                        'temp_limit').get('low'))
                                    T_limit_high = float(each_phase.find(
                                        'temp_limit').get('high'))
                                    try:
                                        each_phase.find('calc_quality').text
                                        calc_quality = str(each_phase.find(
                                            'calc_quality').text)
                                    except:
                                        calc_quality = None
                                    mm = float(each_phase.find(
                                        'molecular_weight').text) / 1e3
                                    coefs = each_phase.find('coefficients')
                                    high_coefs = np.empty((7), dtype='d')
                                    low_coefs = np.empty((7), dtype='d')
                                    range_1000_to_Tmax = coefs.find(
                                        'range_1000_to_Tmax').findall('coef')
                                    for (index, a_term) in enumerate(
                                            range_1000_to_Tmax):
                                        high_coefs[index] = a_term.text
                                    range_Tmin_to_1000 = coefs.find(
                                        'range_Tmin_to_1000').findall('coef')
                                    for (index, a_term) in enumerate(
                                            range_Tmin_to_1000):
                                        low_coefs[index] = a_term.text
                                    h_formation = float(
                                        coefs.find('hf298_div_r').text) * _R
                                    return Compound(
                                        cas,
                                        description,
                                        reference,
                                        formula,
                                        elements,
                                        aggr_state,
                                        T_limit_low,
                                        T_limit_high,
                                        calc_quality,
                                        mm,
                                        low_coefs,
                                        high_coefs,
                                        h_formation)
                    except:
                        pass
            except:
                pass
        return None


# inherits from Elementdb so there is no need to slow down reading
# burcat.xml all the time
class Reaction(Database):
    """Models a simple reaction. Example:\n
    1 N2 + 3 H2 <-> 2 NH3\n
    reaction1 = burcat.Reaction(None, 600, ['n2 ref element',
    'h2 ref element'], ['nh3 anharmonic'], [1, 3], [2])"""
    def __init__(self, p, T, reagents, products, rcoefs=None, pcoefs=None):
        Elementdb.__init__(self)
        self.T = T
        self.p = p
        self.reagents = []
        self.products = []
        self._rcoefs = [abs(z) for z in rcoefs]
        self._pcoefs = [abs(z) for z in pcoefs]

        # error checking
        elements_in_reagents = set()
        for compound in reagents:
            if isinstance(compound, Element):
                self.reagents.append(compound)
                elements_in_reagents.add(tuple(compound._elements))
            elif isinstance(compound, str):
                self.reagents.append(self.set_compound(compound))
                elements_in_reagents.add(tuple(
                    self.set_compound(compound)._elements))
        elements_in_products = set()
        for compound in products:
            if isinstance(compound, Element):
                self.products.append(compound)
                elements_in_products.add(tuple(compound._elements))
            elif isinstance(compound, str):
                self.products.append(self.set_compound(compound))
                elements_in_products.add(tuple(
                    self.set_compound(compound)._elements))

        if set([x[0] for comp in elements_in_products for x in comp]) !=      \
           set([x[0] for comp in elements_in_reagents for x in comp]):
            raise Exception('Cannot balance equation with different'
                            'elements in reagents and products')

        self.deltah = self._delta_enthalpy()
        self.deltag = self._delta_gibbs_energy()
        self.equilibrium_constant = self.equilibrium_constant()

    def _delta_enthalpy(self, T=None):
        """Reaction deltaH in J/mol"""
        if T is not None:
            self.T = T
        delta_h = 0
        for (coefficient, compound) in zip(self._rcoefs, self.reagents):
            delta_h = delta_h - coefficient * compound.enthalpy_engineering(
                self.T)
        for (coefficient, compound) in zip(self._pcoefs, self.products):
            delta_h = delta_h + coefficient * compound.enthalpy_engineering(
                self.T)
        return delta_h

    def _delta_gibbs_energy(self, T=None):
        """Reaction deltaG in J/mol"""
        if T is not None:
            self.T = T
        deltag = 0
        for (coef, comp) in zip(self._rcoefs, self.reagents):
            deltag = deltag - coef * comp.gibbs_energy(self.T)
        for (coef, comp) in zip(self._pcoefs, self.products):
            deltag = deltag + coef * comp.gibbs_energy(self.T)
        return deltag

    def equilibrium_constant(self, T=None):
        """The equilibrium constant: K_eq = exp( - deltaG / (R * T))"""
        if T is not None:
            self.T = T
        return np.exp(-1 * self.deltag / (_R * self.T))

    def __repr__(self):
        r = ''
        for (reag, coef) in zip(self.reagents, self._rcoefs):
            r = r + '+' + str(coef) + ' ' + reag.inp_name + ' '
        r = r + ' -> '
        for (reag, coef) in zip(self.products, self._pcoefs):
            r = r + '+' + str(coef) + ' ' + reag.inp_name + ' '
        return """<reaction> {0}""".format(r)
