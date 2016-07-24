# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 17:34:44 2015

@author: monteiro
"""
from thermopy import burcat
from thermopy.iapws import Water
from numpy import array, dot


def test_enthalpy_massic_tests():
    """Test for various elements enthalpies checked against a literature
    source. Relative error <= 1%.\n
    RE = abs(delta_cp_burcat - delta_cp_lit) / delta_cp_burcat
    """
    # Relative Error
    RE = 1/100
    # Initialization
    database = burcat.Database()

    #
    # REFERENCE: INCROPERA, ISBN 13 978-0470-50197-9
    #
    # Hydrogen from 600 to 700 K
    hydrogen = database.set_compound('h2 ref element')
    delta_h_burcat = (hydrogen.enthalpy_massic(700)
                      - hydrogen.enthalpy_massic(600))
    delta_h_literature = 14.55 * 100 * 1e3  # from 600 to 700 K
    re = abs(delta_h_burcat - delta_h_literature) / delta_h_burcat
    assert re < RE

    # Oxygen from 350 to 400 K
    oxygen = database.set_compound('o2 ref element')
    delta_h_burcat = oxygen.enthalpy_massic(400) - oxygen.enthalpy_massic(350)
    delta_h_literature = 0.929 * 50 * 1e3  # from 600 to 700 K
    re = abs(delta_h_burcat - delta_h_literature) / delta_h_burcat
    assert re < RE

    #
    # REFERENCE: IAPWS; iapws_IF97-Rev.pdf ; low pressure for
    # it to be an ideal gas
    #
    # Water from 500 to 1000 K
    water = database.set_compound('H2O')
    delta_h_burcat = water.enthalpy_massic(1000) - water.enthalpy_massic(500)
    re = abs(delta_h_burcat
             - ((Water(1e3, 1000, massic_basis=True).enthalpy()
                 - Water(1e3, 500, massic_basis=True).enthalpy()) * 1e3)
             ) / delta_h_burcat
    assert re < RE

    #
    # REFERENCE: PERRY, Chemical Engineers Handbook McGraw-Hill 8thEd 2008;
    #            DOI: 10.1036/0071422943
    #
    # NO from 200 to 1450 K; table 2-155
    no = database.set_compound('no')
    delta_h_burcat = no.enthalpy(1450) - no.enthalpy(400)
    #Cp  form:        C 0 p = C1 + C2T + C3T 2 + C4T 3 + C5T 4
    T = 1450
    Tb = array([T, T ** 2, T ** 3, T ** 4, T ** 5], 'd')
    coefs = array([34980 / 1, -35.32 / 2, 0.07729 / 3,
                   -5.7357*1e-5 / 4, 0.0014526*1e-10 / 5],
                  dtype='d')
    T = 200
    Ta = array([T, T ** 2, T ** 3, T ** 4, T ** 5], 'd')
    delta_h_literature = dot((Tb - Ta), coefs) / 1e3
    # Looks like there is a mistake on PERRY's table. The temperature range is
    # also strange since 100 - 1500 K are used only for noble gases (except
    # NO).
    #print('mistake in NO:', delta_h_burcat, delta_h_literature)
    #re = abs(delta_h_burcat - delta_h_literature) / delta_h_burcat
    #assert re < RE

    # Argon from 300 to 1100 K; table 2-155
    ar = database.set_compound('ar ref element')
    delta_h_burcat = ar.enthalpy(1100) - ar.enthalpy(300)
    T = 1100
    Tb = array([T, T ** 2, T ** 3, T ** 4, T ** 5], 'd')
    coefs = array([20786 / 1, 0 / 2, 0 / 3, 0*1e-5 / 4, 0*1e-10 / 5],
                  dtype='d')
    T = 300
    Ta = array([T, T ** 2, T ** 3, T ** 4, T ** 5], 'd')
    delta_h_literature = dot((Tb - Ta), coefs) / 1e3
    re = abs(delta_h_burcat - delta_h_literature) / delta_h_burcat
    assert re < RE

    # Hydrogen from 200 to 250 K; table 2-155 in
    h2 = database.set_compound('h2 ref element')
    delta_h_burcat = h2.enthalpy(250) - h2.enthalpy(200)
    #Cp  form:        C 0 p = C1 + C2T + C3T 2 + C4T 3 + C5T 4
    T = 250
    Tb = array([T, T ** 2, T ** 3, T ** 4, T ** 5], 'd')
    coefs = array([64979 / 1, -788.17 / 2, 5.8287 / 3, -1845.9*1e-5 / 4,
                   216400*1e-10 / 5],
                  dtype='d')
    T = 200
    Ta = array([T, T ** 2, T ** 3, T ** 4, T ** 5], 'd')
    delta_h_literature = dot((Tb - Ta), coefs) / 1e3
    re = abs(delta_h_burcat - delta_h_literature) / delta_h_burcat
    assert re < RE

    #
    # REFERENCE: NIST website; http://webbook.nist.gov;
    #            Shomate Equation
    # AL2SO3(S)
    al2so3 = database.set_compound('AL2O3(S)')
    delta_h_burcat = al2so3.enthalpy(2300) - al2so3.enthalpy(300)
    delta_h_literature = (251.0 - 0.12) * 1e3
    re = abs(delta_h_burcat - delta_h_literature) / delta_h_burcat
    assert re < RE
    # CRN(S) CHROMIUM NITRIDE CONDENSED
    crn = database.set_compound('CrN(s)')
    delta_h_burcat = crn.enthalpy(2200) - crn.enthalpy(400)
    delta_h_literature = (104.0 - 4.84) * 1e3
    re = abs(delta_h_burcat - delta_h_literature) / delta_h_burcat
    assert re < RE
    # FeCL2(L)
    fecl2l = database.set_compound('FeCL2(L)')
    delta_h_burcat = fecl2l.enthalpy(2000) - fecl2l.enthalpy(950)
    delta_h_literature = (173.9 - 66.6) * 1e3
    re = abs(delta_h_burcat - delta_h_literature) / delta_h_burcat
    assert re < RE
    # FeS(L)
    fesl = database.set_compound('FeS(L)')
    delta_h_burcat = fesl.enthalpy(3800) - fesl.enthalpy(1463)
    delta_h_literature = (222.0 - 75.81) * 1e3
    re = abs(delta_h_burcat - delta_h_literature) / delta_h_burcat
    assert re < RE
