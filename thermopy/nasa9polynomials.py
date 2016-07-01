u"""
Enable the NASA 9 term polynomials to create chemical compound objects.

This module allows the user to use the nasa 9 term polynomials (abbreviated as NASA 9) to model chemical compounds and reactions. Is worth noting that all the output is stripped of any phyisical unit, that is, results are returned as numpy floats. Therefore to end any form of ambiguity we reinstate that all the results are in SI units using a molar base. It is up to the user to beware of any physical unit conversion.

Classes:
    Compound: chemical compound present in the NASA9 term polynomials database.

    CompoundIdealGas: chemical compound as an ideal gas present in the NASA9 term polynomials database. It inherits from Compound.

    Reaction: model of a chemical reaction.

Example:
    >>> from thermopy import nasa9polynomials as nasa9
    >>> db = nasa9.Database()
    >>> caf2 = db.set_compound('caf2')
    >>> print(caf2.inchikey)
    WUKWITHWXAAZEY-UHFFFAOYSA-L
    >>> print(caf2.enthalpy_of_formation)
    -790828.409
    >>> print(caf2.heat_capacity(300))
    51.2707324499
    >>> print(caf2.molecular_weight) # note that the SI unit is mol/kg
    0.0780748064
    >>> water = db.set_compound('h2o(l)')
    >>> print(water.entropy(300))
    69.633703
    >>> print(water.elements)
    [('H', '2'), ('O', '1')]

References:
    [1] Bonnie J. McBride, Michael J. Zehe, and Sanford Gordon. NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species. September 2002.

"""

import os
import re
import xml.etree.ElementTree as ET
import numpy as np
from thermopy.constants import ideal_gas_constant
_R = ideal_gas_constant[0]


class Compound(object):
    u"""
    Chemical compound present in the NASA9 term polynomials database.

    It is usually set by nasa9polyniomials.Database.set_compound(xml_compound).

    It has all the thermodynamic functions listed in [1] as methods which take temperature as their sole argument. Those were expanded to include gibbs_energy that could be defined by the given functions.

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
        reference (str): Reference for the compound. See [1] for details.

    Methods:
        enthalpy: calculates the enthalpy for a Compound object.
        entropy:  calculates the entropy for a Compound object.
        gibbs_energy: calculates the Gibbs energy for a Compound object.
        heat_capacity: calculates the heat capacity for a Compound object.

    Subclasses:
        CompoundIdealGas: chemical compound as an ideal gas present in the NASA9 term polynomials database. It inherits from Compound.

    """

    def __init__(self, xml_compound):
        u"""
        Instantiate a Compound object from xml info.

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
        self.elements = self._get_xml_compounds(xml_compound.find(
            'elements'))
        self.condensed = bool(
            xml_compound.find('condensed').text == 'True')
        self.molecular_weight = float(
            xml_compound.find('molecular_weight').text)
        self.enthalpy_of_formation = float(
            xml_compound.find('hf298.15').text)

    def _get_xml_compounds(self, xml_compound):
        u"""
        Return a list of tuples containing the elements and their proportion in the chemical compound.

        Args:
            xml_compound (list): list of Element Tree objects containing elements.

        Returns:
            list: list of tuples. Tuples are of the form (str, int) where str is the symbol of the element and int is its proportion in the molecule.

        """
        xml_compounds_list = []
        for one_xml_compound in xml_compound:
            xml_compounds_list.append(tuple(*map(
                lambda x, y: (x, int(y)),
                *one_xml_compound.items()[0])))
        return xml_compounds_list

    def _evaluate_temperature_interval(self, T):
        u"""
        Ouput temperature interval to be used by public methods.

        Helper method to ouput temperature interval order.

        Args:
            T (float): temperature.

        Returns:
            int: The order of the temperature range (0th, 1st, 2nd, etc). Some compounds have more than one temperature range with different corrisponding coefficients. Therefore the temperature range has to be specified.

        Example:
            The KI gas has two temperature intervals (thus two sets of coefficients to be used). The ranges are: [200, 1000] and [1000, 6000] as for most gases. Thus requiring a property to be measured at 1100 K the second interval should be used and this method shall return the number 1 (as opposed to zero).

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
        for (i, coef) in enumerate(self._xml_compound.findall(
                'T_range')[self._evaluate_temperature_interval(T)]):
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
        for (i, coef) in enumerate(self._xml_compound.findall(
                'T_range')[self._evaluate_temperature_interval(T)]):
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
        for (i, coef) in enumerate(self._xml_compound.findall(
                'T_range')[self._evaluate_temperature_interval(T)]):
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
    Chemical compound as an ideal gas present in the NASA9 term polynomials database.

    Subclasses Compound. It adds two extra methods which come from the thermodynamics of Ideal Gases.

    Ideal gases have an extended set of functions such as internal energy and constant volume heat capacity. All of these properties are derived from definitions on thermodynamics.

    Methods:
        heat_capacity_constant_v: calculates the heat capacity at a constant volume for a CompoundIdealGas object.
        internal_energy: calculates the internal energy for a CompoundIdealGas object.

    """
    def __init__(self, Element):
        u"""Initialize an ideal gas Compound."""
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

    The preferred method for identifying compounds is via InChIKey.
    Other methods are:
        1. CAS number.
        2. Usual name.
        3. IUPAC name.

    Example:
        >>> from thermopy import nasa9polynomials as nasa9
        >>> db = nasa9.Database()

    """

    def __init__(self):
        u"""Initializes the database."""
        xmlPath = os.path.join(
            os.path.dirname(__file__), os.pardir, 'databases',
            'nasa9polynomials.xml')
        self._nasa9 = ET.parse(os.path.abspath(xmlPath))
        self._root = self._nasa9.getroot()

    def _search_database(self, x):
        u"""
        Search the database.

        Args:
            x (str): identifier for compound being searched.

        Returns:
            tuple: (inp file name, iupac name, ET.Element).

        Note:
            The preferred method for identifying compounds is via InChIKey.
            Other methods are:
                1. CAS number.
                2. Usual name.
                3. IUPAC name.

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
                if (x.lower() == specie.find('identification').find(
                            'IUPAC_name').text.lower() or
                        x.lower() == specie.attrib['inp_file_name'].lower()):
                    result_list.append(specie)
            if len(result_list) == 1:
                pass
            else:  # exact match was not sucessfull go to loose match
                result_list = []
                # if not found tries loose match
                for specie in self._root:
                    iupac_name = specie.find('identification').find(
                        'IUPAC_name')
                    augmented_namespace = (specie.attrib['inp_file_name']
                                           + ' ' + specie.find('comment').text
                                           + ' ' + iupac_name.text)
                    if x.lower() in augmented_namespace.lower():
                        result_list.append(specie)
        for specie in result_list:
            return [(y.attrib['inp_file_name'],
                     y.find('identification').find('IUPAC_name').text,
                     y) for y in result_list]

    def list_compound(self, x):
        u"""
        List the compounds for a given input.
        It is intended to be used in interactive mode.

        Args:
            x (str): identifier for compound being searched.

        Returns:
            list: list of tuples containing (str, str) being the 'inp name' and the IUPAC name respectively.

        Note:
            The preferred method for identifying compounds is via InChIKey.
            Other methods are:
                1. CAS number.
                2. Usual name.
                3. IUPAC name.

        Example:
            >>> from thermopy import nasa9polynomials as nasa9
            >>> db = nasa9.Database()
            >>> for i in db.list_compound('fes'):
            ...     print(i)
            ('FeS(a)', 'sulfanylideneiron')
            ('FeS(b)', 'sulfanylideneiron')
            ('FeS(c)', 'sulfanylideneiron')
            ('FeS(L)', 'sulfanylideneiron')
            ('FeSO4(cr)', 'iron(2+);sulfate')
            ('FeS2(cr)', 'N/A')

        """
        result_list = []
        for i in self._search_database(x):
            result_list.append((i[0], i[1]))
        return result_list

    def set_compound(self, x):
        u"""
        Set the compound if there is one entry specified on the database.

        It is important to notice that due to the nature of the work of this database, compounds are gases unless explicitly stated otherwise.

        Args:
            x (str): identifier for compound being searched.

        Returns:
            Compound??fix to thermopy.nasa9polynomials.Compound: returns a Compound objected if the phase is condensed. Returns a CompoundIdealGas otherwise.

        Note:
            The preferred method for identifying compounds is via InChIKey.
            Other methods are:
                1. CAS number.
                2. Usual name.
                3. IUPAC name.

        Example:
            >>> # Someone is looking for the element gallium but is not certain how to instantiate it. One would first list the compounds:
            >>> from thermopy import nasa9polynomials as nasa9
            >>> db = nasa9.Database()
            >>> for i in db.list_compound('gallium'):
            ...     print(i)
            ...
            ('Ga', 'gallium')
            ('Ga+', 'gallium')
            ('GaBr', 'bromogallium')
            ('GaBr2', 'dibromogallium')
            ('GaCl', 'chlorogallium')
            ('GaCl2', 'gallium;dichloride')
            ('GaF2', 'difluorogallium')
            ('GaI2', 'diiodogallium')
            ('GaO', 'oxogallium')
            ('GaOH', 'gallium;hydroxide')
            ('Ga2Cl4', 'gallium;gallium;tetrachloride')
            ('Ga2O', 'gallium;oxygen(2-)')
            ('Ga(cr)', 'gallium')
            ('Ga(L)', 'gallium')
            ('Ga2O3(cr)', 'digallium;oxygen(2-)')
            ('Ga2O3(L)', 'digallium;oxygen(2-)')
            >>> # Finally defining his compound:
            >>> gallium = db.set_compound('Ga')
            >>> gallium.inp_name
            'Ga'

        """
        result = self._search_database(x)
        if len(result) != 1:  # could not set component: give error messages
            if len(result) == 0:
                raise Exception('No compound found.')
            else:
                raise Exception('The compound \'' + str(x) + '\' you are '
                                'trying to set is not unique: ' + result[0][0],
                                result[1][0])
        if result[0][2].find('condensed') == 'True':
            return Compound(result[0][2])
        else:  # if it is an ideal gas
            return CompoundIdealGas(result[0][2])


class Reaction(object):
    u"""
    Model of a chemical reaction.

    """
    def __init__(self, T, reactants, products,
                 reactants_coefficients, product_coefficients):
        u"""
        Initializes a Reaction object.

        Args:
            T (float): temperature of the reaction.
            reactants (tuple): tuple of the reactants as Compounds.
            products (tuple): tuple of the products as Compounds.
            reactants_coefficients (tuple): tuple of the reactants coefficients.
            product_coefficients (tuple): tuple of the products coefficients.

        Example:
            >>> # For the reaction 2 Na + 2 H2O -> 2 NaOH + H2 one would do
            >>> from thermopy import nasa9polynomials as nasa9
            >>> db = nasa9.Database()
            >>> na = db.set_compound('na(cr)')  # being careful to initilize solid compounds
            >>> water = db.set_compound('h2o(l)') # being careful to initilize liquid compounds
            >>> sodium_hydroxide = db.set_compound('naoh(a)') # being careful to initilize solid compounds
            >>> hydrogen = db.set_compound('h2')
            >>> reacts = (na, water)
            >>> prods = (sodium_hydroxide, hydrogen)
            >>> reacts_coefs = (2, 2)
            >>> prods_coefs = (2, 1)
            >>> reaction1 = nasa9.Reaction(300, reacts, prods, reacts_coefs, prods_coefs)
            >>> print(reaction1)
            <reaction> +2 Na(cr) +2 H2O(L)  -> +2 NaOH(a) +1 H2
            >>> print(reaction1.entropy_difference())
            149.097547531
            >>> print(reaction1.enthalpy_reaction())
            -279857.367433

        Note:
            The tuple of the reactants and its coefficients should refer to the same compounds (follow the same order). See example.

        """
        self.T = T
        self._reactants = reactants
        self._products = products
        self._rcoefs = tuple(abs(z) for z in reactants_coefficients)
        self._pcoefs = tuple(abs(z) for z in product_coefficients)
        # error checking
        if (len(self._reactants) != len(self._rcoefs) or
                len(self._products) != len(self._pcoefs)):
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
            deltah = deltah - coefficient * compound.enthalpy(self.T)
        for (coefficient, compound) in zip(self._rcoefs, self._products):
            deltah = deltah + coefficient * compound.enthalpy(self.T)
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
            deltas = deltas - coefficient * compound.entropy(self.T)
        for (coefficient, compound) in zip(self._rcoefs, self._products):
            deltas = deltas + coefficient * compound.entropy(self.T)
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
            deltag = deltag - coefficient * compound.gibbs_energy(self.T)
        for (coefficient, compound) in zip(self._rcoefs, self._products):
            deltag = deltag + coefficient * compound.gibbs_energy(self.T)
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
        return np.exp(-1 * self.gibbs_energy_difference(self.T) / (
            _R * self.T))

    def __repr__(self):
        u"""Define how a reaction should be print."""
        r = ''
        for (reag, coef) in zip(self._reactants, self._rcoefs):
            r = r + '+' + str(coef) + ' ' + reag.inp_name + ' '
        r = r + ' -> '
        for (reag, coef) in zip(self._products, self._pcoefs):
            r = r + '+' + str(coef) + ' ' + reag.inp_name + ' '
        return """<reaction> {0}""".format(r)
