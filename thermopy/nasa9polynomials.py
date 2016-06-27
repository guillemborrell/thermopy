# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 10:00:28 2015

@author: monteiro

Description:

Reads the nasa9polynomials database and returns objects. Input should be of the
type:
benzene = nasa9.get_compound('name_to_search')
name_to_search prefereably should be an InChIKey identifier.
Other possibilities are canonical smiles, cas_number or a simple name (C6H6).
Note that the nasa9polynomials are functions of temperature only.

Then usage should be simple:
benzene.enthalpy(600) - benzene.enthalpy(500)
"""
import re
import xml.etree.ElementTree as ET
import numpy as np
from thermopy.constants import ideal_gas_constant
import os
_R = ideal_gas_constant[0]


class Compound(object):
    
    """
    Chemical compound.
    
    It is usuall set by Database.set_compound(xml_compound).
    
    Attributes:
        canonical_smiles (str): Canonical SMILES (Simplified molecular-input line-entry system) of the compound.
        cas_number (str): CAS (Chemical Abstract Service) number of the compound.
        comment (str): Comment found in the xml database. Usually references.
        condensed (bool): True if the compound is condensed, False if not.
        xml_compounds (list): List containing tuples of two entries. The first is the xml_compound and the second is the proportion of the xml_compound in the molecule. Both values are strings.
        enthalpy_of_formation (float): Enthalpy of formation of the compound.
        inchikey (str): InChI (International Chemical Identifier) key for the compound.
        inp_name (str): ?
        iupac_name (str): IUPAC (International Union of Pure and Applied Chemistry) name of the compound.
        molecular_weight (float): Molecular weight of the compound.
        reference (str): Reference for the compound. See ? for details.
    
    Methods:
        ? one line methods
        
    """






    """The compound set by Database.set_compound(xml_compound).
    Has the functionspresent in the nasa report and some
    other which can be used with the given information.\n
    Nasa 9 functions:\n
    - heat_capacity(T)\n
    - enthalpy(T)
    - entropy(T)
    Extended functions:\n
    - gibbs_energy(T) = H - T * S"""

    def __init__(self, xml_compound):
        """
        
        
        """
        self._xml_compound = xml_compound
        self.inp_name = xml_compound.attrib['inp_file_name']
        self.inchikey = xml_compound.find('identification').find('InChIKey').text
        self.canonical_smiles =                                               \
            xml_compound.find('identification').find('canonical_smiles').text
        self.cas_number =                                                     \
            xml_compound.find('identification').find('cas_number').text
        self.iupac_name =                                                     \
            xml_compound.find('identification').find('IUPAC_name').text
        self.comment = xml_compound.find('comment').text
        self.reference = xml_compound.find('reference').text
        self.xml_compounds = self._get_xml_compounds(xml_compound.find('xml_compounds'))
        self.condensed = bool(xml_compound.find('condensed').text == 'True')
        self.molecular_weight = float(xml_compound.find('molecular_weight').text)
        self.enthalpy_of_formation = float(xml_compound.find('hf298.15').text)

    def _get_xml_compounds(self, xml_compound):
        """
        Helper method that returns a list of tuples containing the xml_compounds and their proportion in the chemical compound.
        
        Args:
            xml_compound
        
        """
        xml_compounds_list = []
        for one_xml_compound in xml_compound:
            xml_compounds_list.append(one_xml_compound.items()[0])
        return(xml_compounds_list)

    def _evaluate_temperature_interval(self, T):
        """Helper method to ouput temperature interval"""
        for (i, Trange) in enumerate(self._xml_compound.findall('T_range')):
            if float(Trange.attrib['Tlow']) <= T <=                           \
               float(Trange.attrib['Thigh']):
                return i
        raise Exception('Temperature out of range for ' + self.iupac_name +
                        '/' + self.inp_name +
                        ' with ' + str(T) + ' K')

    def heat_capacity(self, T):
        """Molar heat capacity at constant pressure at
        temperature T for standard state.\nUnits: J/mol K"""
        coefficients = np.empty(9, dtype=np.float32)
        for (i, coef) in enumerate(self._xml_compound.findall('T_range')[self._evaluate_temperature_interval(T)]):
            coefficients[i] = np.array(coef.text, dtype=np.float32)
        exponents = np.array([-2, -1, 0, 1, 2, 3, 4], dtype=np.signedinteger)
        return np.sum(
               np.multiply(
               np.power(T, exponents, dtype=np.float32),
               coefficients[0:7]), dtype=np.float32) * _R

    def enthalpy(self, T):
        """Molar enthalpy at constant pressure at temperature
        T for standard state.\nUnits: J/mol"""
        coefficients = np.empty(9, dtype=np.float32)
        for (i, coef) in enumerate(self._xml_compound.findall('T_range')[self._evaluate_temperature_interval(T)]):
            coefficients[i] = np.array(coef.text, dtype=np.float32)
        exponents = np.array([-2, -1, 0, 1, 2, 3, 4, -1],
                             dtype=np.signedinteger)
        other_factors = np.array([-1, np.log(T), 1, 0.5, 1/3, 0.25, 0.2, 1],
                                 dtype=np.float32)
        return np.sum(
               np.multiply(
               np.multiply(
               np.power(T, exponents, dtype=np.float32),
               coefficients[0:8]),
               other_factors), dtype=np.float32) * _R * T

    def entropy(self, T):
        """Molar entropy at constant pressure at temperature
        T for standard state.\nUnits: J/mol K"""
        coefficients = np.empty(9, dtype=np.float32)
        for (i, coef) in enumerate(self._xml_compound.findall('T_range')[self._evaluate_temperature_interval(T)]):
            coefficients[i] = np.array(coef.text, dtype=np.float32)
        exponents = np.array([-2, -1, 0, 1, 2, 3, 4, 0, 0],
                             dtype=np.signedinteger)
        other_factors = np.array([-0.5, -1, np.log(T), 1, 0.5, 1/3, 0.25,
                                  0, 1],
                                 dtype=np.float32)
        return np.sum(
               np.multiply(
               np.multiply(
               np.power(T, exponents, dtype=np.float32),
               coefficients[:]),
               other_factors), dtype=np.float32) * _R

    def gibbs_energy(self, T):
        """Molar Gibbs energy at constant pressure at temperature
        T for standard state.\nUnits: J/mol"""
        return self.enthalpy(T) - T * self.entropy(T)

    def __str__(self):
        if self.iupac_name != 'N/A':
            return str(self.iupac_name) + ': ' + self.inp_name
        elif self.canonical_smiles != 'N/A':
            return str(self.canonical_smiles)
        else:
            return str(self.inp_name)


class Compound_ideal_gas(Compound):
    """Subclass of Compounds which have an extended set of functions
    such as internal energy and constant volume heat capacity."""
    def __init__(self, Element):
        Compound.__init__(self, Element)

    def heat_capacity_constant_v(self, T):
            """This function is not explicitly given by the nasa 9
            polynomials. From thermodynamic relations one can infer:\n
            C_v = C_p - R\nMolar heat capacity at constant volune
            at temperature T for standard state.\nUnits: J/mol K"""
            return self.heat_capacity(T) - _R

    def internal_energy(self, T):
        """This function is not explicitly given by the nasa 9
        polynomials. From thermodynamic relations one can
        infer for an ideal gas:\n
        U = H - R * T\nUnits: J/mol"""
        return self.enthalpy(T) - _R * T

class Database(object):
    """
    Nasa 9 polynomials database (see NASA/TP—2002-211556).
    The preferred method for identifying compounds is via InChIKey.
    Other methods are: cas number, usual name and iupac name.

    Example:
    import nasa9polynomials as nasa9

    db = nasa9.Database()

    benzene = db.set_compound('benzene')

    oxygen = db.set_compound('MYMOFIZGZYHOMD-UHFFFAOYSA-N')

    benzene.condensed
    False

    benzene.enthalpy_of_formation
    82880.0

    oxygen.heat_capacity(1000)
    34.882379221984863
    """

    def __init__(self):
        xmlPath = os.path.join(os.path.dirname(__file__),
            os.pardir, 'databases', 'nasa9polynomials.xml')
        self._nasa9 = ET.parse(os.path.abspath(xmlPath))
        self._root = self._nasa9.getroot()

    def _search_database(self, x):
        """Helper function to search database.
        Retuns a list of tuples with 3 terms:\n
        (inp file name,\tiupac name,\tET.Element)"""
        result_list = []
        inchikey_re = re.compile('[A-Z]{14}\-[A-Z]{10}\-[A-Z]')
        cas_re = re.compile('[0-9]{2,7}\-[0-9][0-9]\-[0-9]')
        # InChIKey search
        if re.match(inchikey_re, x):  # is an inchikey
            for specie in self._root:
                identification = specie.find('identification')
                inchikey = identification.find('InChIKey')
                if x == inchikey.text:
                    result_list.append(specie)
        # CAS search
        elif re.match(cas_re, x):
            for specie in self._root:
                identification = specie.find('identification')
                cas_number = identification.find('cas_number')
                if x == cas_number.text:
                    result_list.append(specie)
        # usual name search
        else:
            # tries exact match first
            for specie in self._root:
                iupac_name = specie.find('identification').find('IUPAC_name')
                if x.lower() ==                                               \
                   specie.find('identification').find('IUPAC_name').text.lower() or \
                   x.lower() == specie.attrib['inp_file_name'].lower():
                    result_list.append(specie)
            if len(result_list) == 1:
                pass
            else:  # exact match was not sucessfull go to loose match
                result_list = []
                # if not found tries loose match
                for specie in self._root:
                    iupac_name = specie.find('identification').find('IUPAC_name')
                    augmented_namespace = specie.attrib['inp_file_name'] +    \
                                          ' ' +                               \
                                          specie.find('comment').text + ' ' + \
                                          iupac_name.text
                    if x.lower() in augmented_namespace.lower():
                        result_list.append(specie)
        for specie in result_list:
            return [(y.attrib['inp_file_name'],
                     y.find('identification').find('IUPAC_name').text,
                     y) for y in result_list]

    def list_compound(self, x):
        """For interactive searching of compounds.
        Input can be a cas string, InChIKey or usual name.\n
        Preference is given to exact searches. So for example:\n\n
        db.list_compound('h2o')  returns
        [('H2O', 'oxidane')]\n
        db.list_compound('h2o(')  returns
        [('H2O(cr)', 'oxidane'), ('H2O(L)', 'oxidane')]\n
        On the other hand:\n
        db.list_compound('iron') returns 26 entries all
        of the iupac names containing 'iron'
        """
        result_list = []
        for i in self._search_database(x):
            result_list.append((i[0], i[1]))
        return result_list

    def set_compound(self, x):
        """Sets the compound if there is one entry specified on the database.
        For example looking for liquid water one would specify 'h2o(l)' and for
        steam 'h2o'. Suppose you want FeS(cr) then it would be better to search
        first:\n\n
        db.list_compound('fes')
        [('FeS(a)', 'sulfanylideneiron'),
         ('FeS(b)', 'sulfanylideneiron'),
         ('FeS(c)', 'sulfanylideneiron'),
         ('FeS(L)', 'sulfanylideneiron'),
         ('FeSO4(cr)', 'iron(2+);sulfate'),
         ('FeS2(cr)', 'N/A')]

        fes = db.set_compound('FeS(b)')
        """
        result = self._search_database(x)
        if len(result) != 1:  # could not set component: give error messages
            if len(result) == 0:
                raise Exception('No compound found.')
            else:
                raise Exception('The compound \'' + str(x) +
                                '\' you are trying to set is not unique: ' +
                                result[0][0], result[1][0])
        if result[0][2].find('condensed') == 'True':
            return Compound(result[0][2])
        else:  # if it is an ideal gas
            return Compound_ideal_gas(result[0][2])


class Reaction(object):
    """Models a reaction from the nasa9polynomials.
    Available functions are:\n
    - enthalpy_difference(T)
    - entropy_difference(T)
    - gibbs_energy_difference(T)
    - equilibrium_constant(T)
    """
    def __init__(self, T, reactants, products,
                 reactants_coefficients, product_coefficients):
        self.T = T
        self._reactants = reactants
        self._products = products
        self._rcoefs = [abs(z) for z in reactants_coefficients]
        self._pcoefs = [abs(z) for z in product_coefficients]
        # error checking
        if len(self._reactants) != len(self._rcoefs) or                       \
           len(self._products) != len(self._pcoefs):
            raise Exception('Number of reactants or products is different'
                            'from the number of coefficients given')

    def enthalpy_difference(self, T=None):
        """Enthalpy of the reaction.\nUnits: J/mol"""
        if T is not None:
            self.T = T
        deltah = 0
        for (coefficient, compound) in zip(self._rcoefs, self._reactants):
            deltah = deltah - coefficient *                               \
                compound.enthalpy(self.T)
        for (coefficient, compound) in zip(self._rcoefs, self._products):
            deltah = deltah + coefficient *                               \
                compound.enthalpy(self.T)
        return deltah

    def entropy_difference(self, T=None):
        """Enthalpy of the reaction.\nUnits: J/mol K"""
        if T is not None:
            self.T = T
        deltas = 0
        for (coefficient, compound) in zip(self._rcoefs, self._reactants):
            deltas = deltas - coefficient *                               \
                compound.entropy(self.T)
        for (coefficient, compound) in zip(self._rcoefs, self._products):
            deltas = deltas + coefficient *                               \
                compound.entropy(self.T)
        return deltas

    def gibbs_energy_difference(self, T=None):
        """Gibbs energy of the reaction.\nUnits: J/mol"""
        if T is not None:
            self.T = T
        deltag = 0
        for (coefficient, compound) in zip(self._rcoefs, self._reactants):
            deltag = deltag - coefficient *                               \
                compound.gibbs_energy(self.T)
        for (coefficient, compound) in zip(self._rcoefs, self._products):
            deltag = deltag + coefficient *                               \
                compound.gibbs_energy(self.T)
        return deltag

    def equilibrium_constant(self, T=None):
        """Enthalpy of the reaction:\n
        K = exp(- deltaG / (R T))
        \nUnits: (dimensionless)"""
        if T is not None:
            self.T = T
        return np.exp(-1 * self.gibbs_energy_difference(self.T)
                      / (_R * self.T))

    def __repr__(self):
        r = ''
        for (reag, coef) in zip(self._reactants, self._rcoefs):
            r = r + '+' + str(coef) + ' ' + reag.inp_name + ' '
        r = r + ' -> '
        for (reag, coef) in zip(self._products, self._pcoefs):
            r = r + '+' + str(coef) + ' ' + reag.inp_name + ' '
        return """<reaction> {0}""".format(r)
