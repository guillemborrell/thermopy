"""
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
    u"""
    Chemical compound present in the NASA 9 term polynomials database.

    It is usually set by nasa9polyniomials.Database.set_compound(xml_compound).

    It has all the thermodynamic functions listed in ?ref? as methods which take temperature as their sole argument. Those were expanded to include gibbs_energy that could be defined by the given functions.

    Attributes:
        canonical_smiles (str): Canonical SMILES (Simplified molecular-input line-entry system) of the compound.
        cas_number (str): CAS (Chemical Abstract Service) number of the compound.
        comment (str): Comment found in the xml database. Usually references.
        condensed (bool): True if the compound is condensed, False if not.
        xml_compounds (list): List containing tuples of two entries. The first is the xml_compound and the second is the proportion of the xml_compound in the molecule. Both values are strings.
        enthalpy_of_formation (float): Enthalpy of formation of the compound.
        inchikey (str): InChI (International Chemical Identifier) key for the compound.
        inp_name (str): Name as per the original 'inp' file.
        iupac_name (str): IUPAC (International Union of Pure and Applied Chemistry) name of the compound.
        molecular_weight (float): Molecular weight of the compound.
        reference (str): Reference for the compound. See ? for details.

    """

    def __init__(self, xml_compound):
        u"""
        Instantiate a Compound object from xlm info.

        Args:
            xml_compound: xml tree containing the relevant fields to characterize the attributes and boundaries of temperature for which calculations are valid.

        """
        self._xml_compound = xml_compound
        self.inp_name = xml_compound.attrib['inp_file_name']
        self.inchikey = xml_compound.find('identification').find(
            'InChIKey').text
        self.canonical_smiles = xml_compound.find(
            'identification').find('canonical_smiles').text
        self.cas_number = xml_compound.find('identification').find(
            'cas_number').text
        self.iupac_name = xml_compound.find('identification').find(
            'IUPAC_name').text
        self.comment = xml_compound.find('comment').text
        self.reference = xml_compound.find('reference').text
        self.xml_compounds = self._get_xml_compounds(xml_compound.find(
            'xml_compounds'))
        self.condensed = bool(
            xml_compound.find('condensed').text == 'True')
        self.molecular_weight = float(
            xml_compound.find('molecular_weight').text)
        self.enthalpy_of_formation = float(
            xml_compound.find('hf298.15').text)

    def _get_xml_compounds(self, xml_compound):
        u"""
        Helper method that returns a list of tuples containing the elements and their proportion in the chemical compound.

        Args:
            xml_compound

        """
        xml_compounds_list = []
        for one_xml_compound in xml_compound:
            xml_compounds_list.append(one_xml_compound.items()[0])
        return(xml_compounds_list)

    def _evaluate_temperature_interval(self, T):
        u"""
        Helper method to ouput temperature interval order.

        Args:
            T (float): temperature.

        Returns:
            int: The order of the temperature range (0th, 1st, 2nd, etc). Some compounds have more than one temperature range with different corrisponding coefficients. Therefore the temperature range has to be specified.

        """
        for (i, Trange) in enumerate(self._xml_compound.findall('T_range')):
            if (float(Trange.attrib['Tlow']) <= T <= float(
                    Trange.attrib['Thigh'])):
                return i
        raise Exception('Temperature out of range for '
                        + self.iupac_name + '/' + self.inp_name
                        + ' with ' + str(T) + ' K')

    def heat_capacity(self, T):
        u"""
        Calculate molar heat capacity at constant pressure for standard state.

        Args:
            T (float): temperature.

        Returns:
            numpy_float: The molar heat capacity for the compound for a given temperature.

        """
        coefficients = np.empty(9, dtype=np.float32)
        for (i, coef) in enumerate(self._xml_compound.findall('T_range')[self._evaluate_temperature_interval(T)]):
            coefficients[i] = np.array(coef.text, dtype=np.float32)
        exponents = np.array([-2, -1, 0, 1, 2, 3, 4], dtype=np.signedinteger)
        return np.sum(
               np.multiply(
               np.power(T, exponents, dtype=np.float32),
               coefficients[0:7]), dtype=np.float32) * _R

    def enthalpy(self, T):
        u"""
        Calculate molar enthalpy at constant pressure for the compound for a given temperature.
        
        Args:
            T (float): temperature.

        Returns:
            numpy_float: The molar enthalpy for the compound for a given temperature.
            
        """
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
        u"""
        Calculate molar entropy at constant pressure for the compound for a given temperature.
        
        Args:
            T (float): temperature.

        Returns:
            numpy_float: The molar entropy for the compound for a given temperature.
            
        """
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
        u"""
        Calculate molar Gibbs energy at constant pressure for the compound for a given temperature.
        
        Args:
            T (float): temperature.

        Returns:
            numpy_float: The molar entropy for the compound for a given temperature.
            
        """
        return self.enthalpy(T) - T * self.entropy(T)

    def __str__(self):
        if self.iupac_name != 'N/A':
            return str(self.iupac_name) + ': ' + self.inp_name
        elif self.canonical_smiles != 'N/A':
            return str(self.canonical_smiles)
        else:
            return str(self.inp_name)


class CompoundIdealGas(Compound):
    u"""
    Chemical compound as an ideal gas present in the NASA 9 term polynomials database.
    
    Subclasses Compound.
    
    
    Ideal gases have an extended set of functions such as internal energy and constant volume heat capacity. All of these properties are derived from definitions on thermodynamics.
    
    """
    def __init__(self, Element):
        Compound.__init__(self, Element)

    def heat_capacity_constant_v(self, T):
        u"""
        Calculate molar heat capacity at constant volume for standard state.

        Args:
            T (float): temperature.

        Returns:
            numpy_float: The molar heat capacity for the compound for a given temperature.

        """
        return self.heat_capacity(T) - _R

    def internal_energy(self, T):
        u"""
        Calculate molar internal energy at constant pressure for standard state.

        Args:
            T (float): temperature.

        Returns:
            numpy_float: The molar internal energy for the compound for a given temperature.

        """
        return self.enthalpy(T) - _R * T

class Database(object):
    u"""
    Nasa 9 term polynomials database (see NASA/TP—2002-211556).
    
    The preferred method for identifying compounds is via InChIKey. Other methods are: cas number, usual name and iupac name.

    Example:
    
    """

    def __init__(self):
        u"""Initializes the database."""
        xmlPath = os.path.join(os.path.dirname(__file__),
            os.pardir, 'databases', 'nasa9polynomials.xml')
        self._nasa9 = ET.parse(os.path.abspath(xmlPath))
        self._root = self._nasa9.getroot()

    def _search_database(self, x):
        u"""
        Helper function to search database.
        
        Args:
            x ():
                
        Returns:
            tuple: (inp file name, iupac name, ET.Element)
            
        """
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
        u"""
        List the compounds find gor a given input.
        It is intended to be used in interactive mode.
        
        
        Input can be a cas string, InChIKey or usual name. Preference is given to exact searches.
        For example:
        ??
        db.list_compound('h2o')  returns
        [('H2O', 'oxidane')]\n
        db.list_compound('h2o(')  returns
        [('H2O(cr)', 'oxidane'), ('H2O(L)', 'oxidane')]\n
        On the other hand:\n
        db.list_compound('iron') returns 26 entries all
        of the iupac names containing 'iron'
        ??
        """
        result_list = []
        for i in self._search_database(x):
            result_list.append((i[0], i[1]))
        return result_list

    def set_compound(self, x):
        u"""
        Set the compound if there is one entry specified on the database.
        
        It is important to notice that due to the nature of the work of this database, compounds are gases unless explicitly stated otherwise.
        
        Example:
            ??
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
            ??
        
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
    u"""
    Model of a chemical reaction.
    
    """
    def __init__(self, T, reactants, products,
                 reactants_coefficients, product_coefficients):
        u""""""
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

    def enthalpy_reaction(self, T=None):
        u"""
        Calculate the enthalpy of the reaction at the standard state.

        Args:
            T (float): temperature.

        Returns:
            numpy_float: The enthalpy of the reaction for a given temperature.

        """
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
        u"""
        Calculate the entropy of the reaction at the standard state.

        Args:
            T (float): temperature.

        Returns:
            numpy_float: The entropy of the reaction for a given temperature.

        """
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
        u"""
        Calculate the Gibbs energy of the reaction at the standard state.

        Args:
            T (float): temperature.

        Returns:
            numpy_float: The Gibbs energy of the reaction for a given temperature.

        """
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
        u"""
        Calculate the equilibrium constant of the reaction at the standard state.
        
        Definition: K = exp(- deltaG / (R T))

        Args:
            T (float): temperature.

        Returns:
            numpy_float: The Gibbs energy of the reaction for a given temperature.
        
        """
        if T is not None:
            self.T = T
        return np.exp(-1 * self.gibbs_energy_difference(self.T)
                      / (_R * self.T))

    def __repr__(self):
        u"""Define how a reaction should be print."""
        r = ''
        for (reag, coef) in zip(self._reactants, self._rcoefs):
            r = r + '+' + str(coef) + ' ' + reag.inp_name + ' '
        r = r + ' -> '
        for (reag, coef) in zip(self._products, self._pcoefs):
            r = r + '+' + str(coef) + ' ' + reag.inp_name + ' '
        return """<reaction> {0}""".format(r)
