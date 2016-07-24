u"""
Enable the NASA 9 term polynomials to create chemical compound objects.

This module allows the user to use the NASA9 term polynomials (abbreviated as
NASA 9) to model chemical compounds and reactions. It is worth noting that all
the output is stripped of any phyisical unit, that is, results are returned as
numpy floats. *Therefore to end any form of ambiguity we reinstate that all the
results are in SI units using molar basis*. It is up to the user to beware of
any physical unit conversion concerning his problem.

Classes:

    Compound: chemical compound present in the NASA9 term polynomials
    database.

    CompoundIdealGas: chemical compound as an ideal gas present in the NASA9
    term polynomials database. It inherits from Compound.

    Reaction: model of a chemical reaction.

Example:
    >>> import thermopy
    >>> from thermopy import nasa9polynomials as nasa9
    >>> db = nasa9.Database()
    >>> caf2 = db.set_compound('caf2')
    >>> print(caf2.elements)
    [('C', 1), ('F', 2)]
    >>> print(caf2.inchikey)
    WUKWITHWXAAZEY-UHFFFAOYSA-L
    >>> print(caf2.enthalpy_of_formation)
    -790828.409
    >>> print(caf2.heat_capacity(300))
    51.2707324499
    >>> print(caf2.molecular_weight)
    0.0780748064
    >>> water = db.set_compound('h2o(l)')
    >>> print(water.entropy(300))
    69.633703
    >>> print(water.elements)
    [('H', 2), ('O', 1)]

References: [1] Bonnie J. McBride, Michael J. Zehe, and Sanford Gordon. NASA
Glenn Coefficients for Calculating Thermodynamic Properties of Individual
Species. September 2002.

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

    It is usually set by nasa9polyniomials.Database.set_compound('identifier').
    If set this way this method already instantiates either a Compound or a
    CompoundIdealGas based on the phase of the chemical compound specified by
    the identifier (see note below).

    Note: *The default phase for compounds in this database is gas.* Thus if
    instantiating with the name 'H2O' one would get steam. For liquid water and
    ice one would rather look for 'H2O(L)' or 'H2O(cr)'. The same is valid for
    a lot of compounds expected to be in the condensed form (such as NaCl,
    Tungsten, etc).

    It has all the thermodynamic functions listed in [1] as methods which take
    temperature as their sole argument. Those were expanded to include
    gibbs_energy that could be defined by the given functions.

    Attributes:
        canonical_smiles (str): Canonical SMILES
            (Simplified molecular-input line-entry system) of the compound.
        cas_number (str): CAS (Chemical Abstract Service) number of the
            compound.
        comment (str): Comment found in the xml database. Usually references.
        condensed (bool): True if the compound is condensed, False if not.
        xml_compounds (list): List containing tuples of two entries. The first
            is the xml_compound and the second is the proportion of the
            xml_compound in the molecule. Both values are strings.
        enthalpy_of_formation (float): Enthalpy of formation of the compound.
        inchikey (str): InChI (International Chemical Identifier) key
            for the compound.
        inp_name (str): Name as per the original 'inp' file.
        iupac_name (str): IUPAC (International Union of Pure and Applied
            Chemistry) name of the compound.
        molecular_weight (float): Molecular
            weight of the compound.
        reference (str): Reference for the compound. See [1] for details.

    Methods:
        enthalpy:calculates the enthalpy for a Compound object.
        entropy: calculates the entropy for a Compound object.
        gibbs_energy: calculates the Gibbs energy for a Compound object.
        heat_capacity: calculates the heat capacity for a Compound object.

    Subclasses:
        CompoundIdealGas: chemical compound as an ideal gas present in
        the NASA9 term polynomials database. It inherits from Compound.

    Examples:
        >>> import thermopy
        >>> from thermopy import nasa9polynomials as nasa9
        >>> db = nasa9.Database()
        >>> uf6 = db.set_compound('uf6(cr)')
        >>> print(uf6)
        hexafluorouranium: UF6(cr)

    """

    def __init__(self, xml_compound):
        u"""
        Instantiate a Compound object from xml info.

        Arguments:
            xml_compound: xml tree containing the relevant fields to
        characterize the attributes and boundaries of temperature for which
        calculations are valid.

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
        Return a list of tuples containing the elements and their
        proportion in the chemical compound.

        Arguments:
            xml_compound (list): list of Element Tree objects containing
        elements.

        Returns:
            list: list of tuples. Tuples are of the form (str, int) where
        str is the symbol of the element and int is its proportion in the
        molecule.

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

        Arguments:
            T (float): temperature.

        Returns:
            int: The order of the temperature range (0th, 1st, 2nd, etc).
        Some compounds have more than one temperature range with different
        corrisponding coefficients. Therefore the temperature range has to be
        specified.

        Example: The KI gas has two temperature intervals (thus two sets of
        coefficients to be used). The ranges are: [200, 1000] and [1000, 6000]
        as for most gases. Thus requiring a property to be measured at 1100 K
        the second interval should be used and this method shall return the
        number 1 (as opposed to zero).

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
        Calculate molar heat capacity at constant pressure for standard
        state.

        Arguments:
            T (float): temperature.

        Returns:
            numpy_float: The molar heat capacity for the compound for a
        given temperature.

        Examples:
            >>> import thermopy
            >>> from thermopy import nasa9polynomials as nasa9
            >>> db = nasa9.Database()
            >>> libr = db.set_compound('LiBr')
            >>> print(libr, '-', libr.heat_capacity(2700))
            lithium;bromide: LiBr - 39.6593057506

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
        Calculate molar enthalpy at constant pressure for the compound for
        a given temperature.

        Arguments:
            T (float): temperature.

        Returns:
            numpy_float: The molar enthalpy for the compound for a given
            temperature.

        Examples:
            >>> import thermopy
            >>> from thermopy import nasa9polynomials as nasa9
            >>> db = nasa9.Database()
            >>> magnesium_hydroxide = db.set_compound('Mg(OH)2(cr)')
            >>> print(magnesium_hydroxide, '-',
            >>> magnesium_hydroxide.enthalpy(500))  # J/mol
            Mg(OH)2(cr) - -906097.801815

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
        Calculate molar entropy at constant pressure for the compound for
        a given temperature.

        Arguments:
            T (float): temperature.

        Returns:
            numpy_float: The molar entropy for the compound for a given
            temperature.

        Examples:
            >>> import thermopy
            >>> from thermopy import nasa9polynomials as nasa9
            >>> db = nasa9.Database()
            >>> argon = db.set_compound('ar')
            >>> print(argon, '-', argon.entropy(200))
            argon: Ar - 146.546470215

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
                np.multiply(np.power(T, exponents, dtype=np.float32),
                            coefficients[:]),
                other_factors), dtype=np.float32) * _R

    def gibbs_energy(self, T):
        u"""
        Calculate molar Gibbs energy at constant pressure for the compound
        for a given temperature.

        Arguments:
            T (float): temperature.

        Returns:
            numpy_float: The molar entropy for the compound for a given
            temperature.

        Examples:
            >>> import thermopy
            >>> from thermopy import nasa9polynomials as nasa9
            >>> db = nasa9.Database()
            >>> csbr = db.set_compound('csbr(cr)')
            >>> print(csbr, '-', csbr.gibbs_energy(273.15))
            cesium;bromide: CsBr(cr) - -436504.410044

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
    Chemical compound as an ideal gas present in the NASA9 term
    polynomials database.

    It is usually set by nasa9polyniomials.Database.set_compound('identifier').
    If set this way this method already instantiates either a Compound or a
    CompoundIdealGas based on the phase of the chemical compound specified by
    the identifier (see note below).

    Note: *The default phase for compounds in this database is gas.* Thus if
    instantiating with the name 'H2O' one would get steam. For liquid water and
    ice one would rather look for 'H2O(L)' or 'H2O(cr)'. The same is valid for
    a lot of compounds expected to be in the condensed form (such as NaCl,
    Tungsten, etc).

    It adds two extra methods which come from the thermodynamics of Ideal
    Gases.

    Inherits from Compound.

    Methods:
        heat_capacity_constant_v: calculates the heat capacity at a
        constant volume for a CompoundIdealGas object.
        internal_energy: calculates
        the internal energy for a CompoundIdealGas object.

    Examples:
        >>> # Instantiating a Compound whose condensed attributed is False
        >>> # automatically sets it as an Ideal Gas:
        >>> import thermopy
        >>> from thermopy import nasa9polynomials as nasa9
        >>> db = nasa9.Database()
        >>> co2 = db.set_compound('CO2')
        >>> print(co2, type(co2))
        carbon dioxide: CO2 <class
        'thermopy.nasa9polynomials.CompoundIdealGas'>

    """

    def __init__(self, xml_compound):
        u"""
        Initialize an ideal gas Compound.

        Arguments:
            xml_compound: xml tree containing the relevant fields to
        characterize the attributes and boundaries of temperature for which
        calculations are valid.

        """
        Compound.__init__(self, xml_compound)

    def heat_capacity_constant_v(self, T):
        u"""
        Calculate molar heat capacity at constant volume for standard
        state.

        Arguments:
            T (float): temperature.

        Returns:
            numpy_float: The molar heat capacity for the compound for a
            given temperature.

        Examples:
            >>> import thermopy
            >>> from thermopy import nasa9polynomials as nasa9
            >>> db = nasa9.Database()
            >>> xenon = db.set_compound('Xe')
            >>> print(xenon, '-', xenon.heat_capacity(298.15))
            xenon: Xe - 20.78618
            >>> print(xenon, '-', xenon.heat_capacity_constant_v(298.15))
            xenon: Xe - 12.471708
            >>> print('subtracting both hc:',
            ...       xenon.heat_capacity(298.15)
            ...       - xenon.heat_capacity_constant_v(298.15))
            subtracting both hc: 8.314472

        """
        return self.heat_capacity(T) - _R

    def internal_energy(self, T):
        u"""
        Calculate molar internal energy at constant pressure for standard
        state.

        Arguments:
            T (float): temperature.

        Returns:
            numpy_float: The molar internal energy for the compound for a
            given temperature.

        Examples:
            >>> import thermopy
            >>> from thermopy import nasa9polynomials as nasa9
            >>> db = nasa9.Database()
            >>> hcn = db.set_compound('hcn')
            >>> print(hcn, '-', hcn.internal_energy(500))
            formonitrile: HCN - 136809.830971

        """
        return self.enthalpy(T) - _R * T


class Database(object):
    u"""
    Nasa 9 term polynomials database (see NASA/TP—2002-211556).

    The preferred method for identifying compounds is via *usual name* followed
    by an aggregation state if the compound is not a gas. E.g. '(L)', '(cr)',
    '(a)', '(b)' where '(a)' and '(b)' are for allotropes.
    Other methods are:
        1. InChIKey.
        2. CAS number.
        3. IUPAC name.

    Examples:
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

        Arguments:
            x (str): identifier for compound being searched.

        Returns:
            tuple: (inp file name, iupac name, ET.Element).

        Note:
            The preferred method for identifying compounds is via *usual name*
            followed by an aggregation state if the compound is not a gas. E.g.
            '(L)', '(cr)', '(a)', '(b)' where '(a)' and '(b)' are for
            allotropes.
            Other methods are:
                1. InChIKey.
                2. CAS number.
                3. IUPAC name.

        """
        result_list = []
        inchikey_re = re.compile('[A-Z]{14}-[A-Z]{10}-[A-Z]')
        cas_re = re.compile('[0-9]{2,7}-[0-9][0-9]-[0-9]')
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

        Arguments:
            x (str): identifier for compound being searched.

        Returns:
            list: list of tuples containing (str, str) being the 'inp name' and
            the IUPAC name respectively.

        Examples:
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

        Note:
            The preferred method for identifying compounds is via *usual name*
            followed by an aggregation state if the compound is not a gas. E.g.
            '(L)', '(cr)', '(a)', '(b)' where '(a)' and '(b)' are for
            allotropes.
            Other methods are:
                1. InChIKey.
                2. CAS number.
                3. IUPAC name.

        """
        result_list = []
        for i in self._search_database(x):
            result_list.append((i[0], i[1]))
        return result_list

    def set_compound(self, x):
        u"""
        Set the compound if there is one entry specified on the database.

        It is important to notice that due to the nature of the work of this
        database, compounds are gases unless explicitly stated otherwise.

        Arguments:
            x (str): identifier for compound being searched.

        Returns:
            Compound: returns a Compound object if the phase is condensed.
            Returns a CompoundIdealGas otherwise.

        Example:
            >>> # Someone is looking for the element gallium but is not certain
            >>> # how to instantiate it. One would first list the compounds:
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
            >>> gallium = db.set_compound('Ga')
            >>> print(gallium)
            gallium: Ga

        Note:
            The preferred method for identifying compounds is via *usual name*
            followed by an aggregation state if the compound is not a gas. E.g.
            '(L)', '(cr)', '(a)', '(b)' where '(a)' and '(b)' are for
            allotropes.
            Other methods are:
                1. InChIKey.
                2. CAS number.
                3. IUPAC name.

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

    Model of a chemical reaction using the NASA9 Compounds.

    Inherits from object.

    Methods:
        enthalpy_reaction: enthalpy of the reaction.
        entropy_reaction: entropy of the reaction.
        gibbs_energy_reaction: Gibbs energy of the reaction.
        equilibrium_constant: equilibrium constant of the reaction.

    Examples:
        >>> from thermopy import nasa9polynomials as nasa9
        >>> db = nasa9.Database()
        >>> na = db.set_compound('na(cr)')
        >>> # being careful to initilize solid compounds
        >>> water = db.set_compound('h2o(l)')
        >>> # being careful to initilize liquid compounds
        >>> sodium_hydroxide = db.set_compound('naoh(a)')
        >>> # being careful to initilize solid compounds
        >>> hydrogen = db.set_compound('h2')
        >>> reacts = (na, water)
        >>> prods = (sodium_hydroxide, hydrogen)
        >>> reacts_coefs = (2, 2)
        >>> prods_coefs = (2, 1)
        >>> reaction1 = nasa9.Reaction(300, reacts, prods, reacts_coefs,
        >>> prods_coefs)
        >>> print(reaction1)
        <reaction> +2 Na(cr) +2 H2O(L)  -> +2 NaOH(a) +1 H2
        >>> print(reaction1.entropy_reaction())
        149.097547531
        >>> print(reaction1.enthalpy_reaction())
        -279857.367433

    Notes:
        The Reaction class does not check for imbalances of the reaction (yet).

    """

    def __init__(self, T, reactants, products,
                 reactants_coefficients, product_coefficients):
        u"""
        Initializes a Reaction object.

        Arguments:
            T (float): temperature of the reaction.
            reactants (tuple): tuple of the reactants as Compounds.
            products (tuple): tuple of the products as Compounds.
            reactants_coefficients (tuple): tuple of the reactants
                coefficients.
            product_coefficients (tuple): tuple of the products coefficients.

        Examples:
            >>> import thermopy
            >>> from thermopy import nasa9polynomials as nasa9
            >>> db = nasa9.Database()
            >>> co = db.set_compound('carbon monoxide')
            >>> molybdenium_oxide = db.set_compound('MoO2(cr)')  # to set it a
            >>> solid
            >>> co2 = db.set_compound('co2')
            >>> molybdenium = db.set_compound('Mo(cr)')  # to set it a solid
            >>> reagents = (co, molybdenium_oxide)
            >>> products = (co2, molybdenium)
            >>> reactants_stoichometry = (2, 1)
            >>> prodcuts_stoichometry = (2, 1)
            >>> reaction1 = nasa9.Reaction(
            ...     298,
            ...     reagents,
            ...     products,
            ...     reactants_stoichometry,
            ...     prodcuts_stoichometry
            ...     )
            >>> print(reaction1, reaction1.enthalpy_reaction())
            <reaction> +2 CO +1 MoO2(cr)  -> +2 CO2 +1 Mo(cr)  23352.3949968

        Note:
            The tuple of the reactants and its coefficients should refer to the
            same compounds (follow the same order). See example.

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

        Arguments:
            T (float): temperature.

        Returns:
            numpy_float: The enthalpy of the reaction for a given temperature.

        Examples:
            >>> import thermopy
            >>> from thermopy import nasa9polynomials as nasa9
            >>> db = nasa9.Database()
            >>> nitric_acid = db.set_compound('hno3')
            >>> # the liquid phase is not present in the database
            >>> naoh = db.set_compound('naoh(a)')
            >>> sodium_nitrate = db.set_compound('nano3(a)')
            >>> water = db.set_compound('h2o(l)')
            >>> reagents = (nitric_acid, naoh)
            >>> products = (sodium_nitrate, water)
            >>> reactants_stoichometry = (1, 1)
            >>> prodcuts_stoichometry = (1, 1)
            >>> reaction1 = nasa9.Reaction(
            ...     298,
            ...     reagents,
            ...     products,
            ...     reactants_stoichometry,
            ...     prodcuts_stoichometry
            ...     )
            >>> print(reaction1, reaction1.enthalpy_reaction())
            <reaction> +1 HNO3 +1 NaOH(a)  -> +1 NaNO3(a) +1 H2O(L)
            -193773.133358

        """
        if T is not None:
            self.T = T
        deltah = 0
        for (coefficient, compound) in zip(self._rcoefs, self._reactants):
            deltah = deltah - coefficient * compound.enthalpy(self.T)
        for (coefficient, compound) in zip(self._rcoefs, self._products):
            deltah = deltah + coefficient * compound.enthalpy(self.T)
        return deltah

    def entropy_reaction(self, T=None):
        u"""
        Calculate the entropy of the reaction at the standard state.

        Arguments:
            T (float): temperature.

        Returns:
            numpy_float: The entropy of the reaction for a given temperature.

        Examples:
            >>> import thermopy
            >>> from thermopy import nasa9polynomials as nasa9
            >>> db = nasa9.Database()
            >>> nitric_acid = db.set_compound('hno3')
            >>> # the liquid phase is not present in the database
            >>> naoh = db.set_compound('naoh(a)')
            >>> sodium_nitrate = db.set_compound('nano3(a)')
            >>> water = db.set_compound('h2o(l)')
            >>> reagents = (nitric_acid, naoh)
            >>> products = (sodium_nitrate, water)
            >>> reactants_stoichometry = (1, 1)
            >>> prodcuts_stoichometry = (1, 1)
            >>> reaction1 = nasa9.Reaction(
            ...     298,
            ...     reagents,
            ...     products,
            ...     reactants_stoichometry,
            ...     prodcuts_stoichometry
            ...     )
            >>> print(reaction1, reaction1.entropy_reaction())
            <reaction> +1 HNO3 +1 NaOH(a)  -> +1 NaNO3(a) +1 H2O(L)
            -145.797754143

        """
        if T is not None:
            self.T = T
        deltas = 0
        for (coefficient, compound) in zip(self._rcoefs, self._reactants):
            deltas = deltas - coefficient * compound.entropy(self.T)
        for (coefficient, compound) in zip(self._rcoefs, self._products):
            deltas = deltas + coefficient * compound.entropy(self.T)
        return deltas

    def gibbs_energy_reaction(self, T=None):
        u"""
        Calculate the Gibbs energy of the reaction at the standard state.

        Arguments:
            T (float): temperature.

        Returns:
            numpy_float: The Gibbs energy of the reaction for a given
            temperature.

        Examples:
            >>> import thermopy
            >>> from thermopy import nasa9polynomials as nasa9
            >>> db = nasa9.Database()
            >>> pcl5 = db.set_compound('pcl5')
            >>> pcl3 = db.set_compound('pcl3')
            >>> chlorine = db.set_compound('cl2')
            >>> reagents = (pcl5,)
            >>> products = (pcl3, chlorine)
            >>> reactants_stoichometry = (1,)
            >>> prodcuts_stoichometry = (1, 1)
            >>> reaction1 = nasa9.Reaction(
            ...     298,
            ...     reagents,
            ...     products,
            ...     reactants_stoichometry,
            ...     prodcuts_stoichometry
            ...     )
            >>> print(reaction1, reaction1.gibbs_energy_reaction())
            <reaction> +1 PCl5  -> +1 PCl3 +1 Cl2  103038.712535

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
        Calculate the equilibrium constant of the reaction at the standard
        state.

        Definition: K = exp(- deltaG / (R T))

        Arguments:
            T (float): temperature.

        Returns:
            numpy_float: The Gibbs energy of the reaction for a given
            temperature.


        Examples:
            >>> import thermopy
            >>> from thermopy import nasa9polynomials as nasa9
            >>> db = nasa9.Database()
            >>> pcl5 = db.set_compound('pcl5')
            >>> pcl3 = db.set_compound('pcl3')
            >>> chlorine = db.set_compound('cl2')
            >>> reagents = (pcl5,)
            >>> products = (pcl3, chlorine)
            >>> reactants_stoichometry = (1,)
            >>> prodcuts_stoichometry = (1, 1)
            >>> reaction1 = nasa9.Reaction(
            ...     500,
            ...     reagents,
            ...     products,
            ...     reactants_stoichometry,
            ...     prodcuts_stoichometry
            ...     )
            >>> print(reaction1, reaction1.equilibrium_constant())
            <reaction> +1 PCl5  -> +1 PCl3 +1 Cl2  6.39431126134e-13

        """
        if T is not None:
            self.T = T
        return np.exp(-1 * self.gibbs_energy_reaction(self.T) / (
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
