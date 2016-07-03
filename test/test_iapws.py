# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 09:29:56 2015

@author: monteiro
"""
from thermopy.iapws import Water
from thermopy.units import Pressure, Temperature


def test_iapws():
    """
    Tests are given inside the IAPWS document. See references for more details.
    """
    #test Tsat given P
    assert round(Water(
        1e5, 373.15).temperature_saturation(0.1e6), 6) == 372.755919
    assert round(Water(
        1e5, 373.15).temperature_saturation(1e6), 6) == 453.035632
    assert round(Water(
        1e5, 373.15).temperature_saturation(10e6), 6) == 584.149488

    #test Psat given T
    assert round(Water(
        1e5, 373.15).pressure_saturation(300).MPa, 11) == 0.00353658941
    assert round(Water(
        1e5, 373.15).pressure_saturation(500).MPa, 8) == 2.63889776
    assert round(Water(
        1e5, 373.15).pressure_saturation(600).MPa, 7) == 12.3443146

    #test regions
        # arbitrary points
    point_in_region1 = (Pressure(20e6), Temperature(300))
    point_in_region2 = (Pressure(1e5), Temperature(373.15))
    point_in_region3 = (Pressure(40e6), Temperature(700))
    point_in_region4 = (Pressure(1).unit('atm'), Temperature(373.1243))
    point_in_region5 = (Pressure(20e6), Temperature(1500))
    assert Water(*point_in_region1)._is_in_region() == 1
    assert Water(*point_in_region2)._is_in_region() == 2
    assert Water(*point_in_region3)._is_in_region() == 3
    # region 4 does not exist as a region; it is rather the saturation line
    assert Water(*point_in_region5)._is_in_region() == 5
#region 1
    #assert specific volume
    assert round(Water(3e6, 300, massic_basis=True).specific_volume(),
                 11) == 0.00100215168
    assert round(Water(80e6, 300, massic_basis=True).specific_volume(),
                 12) == 0.000971180894
    assert round(Water(3e6, 500, massic_basis=True).specific_volume(),
                 11) == 0.00120241800
#
#    #assert internal energy
    assert round(Water(3e6, 300, massic_basis=True).internal_energy(),
                 6) == 112.324818
    assert round(Water(80e6, 300, massic_basis=True).internal_energy(),
                 6) == 106.448356
    assert round(Water(3e6, 500, massic_basis=True).internal_energy(),
                 6) == 971.934985
#
#    #assert enthropy
    assert round(Water(3e6, 300, massic_basis=True).entropy(),
                 9) == 0.392294792
    assert round(Water(80e6, 300, massic_basis=True).entropy(),
                 9) == 0.368563852
    assert round(Water(3e6, 500, massic_basis=True).entropy(),
                 8) == 2.58041912

    #assert enthalpy
    assert round(Water(3e6, 300, massic_basis=True).enthalpy(),
                 6) == 115.331273
    assert round(Water(80e6, 300, massic_basis=True).enthalpy(),
                 6) == 184.142828
    assert round(Water(3e6, 500, massic_basis=True).enthalpy(),
                 6) == 975.542239

    #assert cp
    assert round(Water(3e6, 300, massic_basis=True).heat_capacity(),
                 8) == 4.17301218
    assert round(Water(80e6, 300, massic_basis=True).heat_capacity(),
                 8) == 4.01008987
    assert round(Water(3e6, 500, massic_basis=True).heat_capacity(),
                 8) == 4.65580682

#    #assert cv
#    assert round(Water(3e6, 300).enthalpy(),6) == 115.331273
#    assert round(Water(80e6, 300).enthalpy(),6) == 184.142828
#    assert round(Water(3e6, 500).enthalpy(),6) == 975.542239
#
    #assert speed of sound
    assert round(Water(3e6, 300, massic_basis=True).speed_of_sound(),
                 5) == 1507.73921
    assert round(Water(80e6, 300, massic_basis=True).speed_of_sound(),
                 5) == 1634.69054
    assert round(Water(3e6, 500, massic_basis=True).speed_of_sound(),
                 5) == 1240.71337

#region 2
    #assert specific volume
    assert round(Water(3500, 300, massic_basis=True).specific_volume(),
                 7) == 39.4913866
    assert round(Water(3500, 700, massic_basis=True).specific_volume(),
                 7) == 92.3015898
    assert round(Water(30e6, 700, massic_basis=True).specific_volume(),
                 11) == 0.00542946619
#
#    #assert internal energy
    assert round(Water(3500, 300, massic_basis=True).internal_energy(),
                 5) == 2411.69160
    assert round(Water(3500, 700, massic_basis=True).internal_energy(),
                 5) == 3012.62819
    assert round(Water(30e6, 700, massic_basis=True).internal_energy(),
                 5) == 2468.61076
#
#    #assert enthropy
    assert round(Water(3500, 300, massic_basis=True).entropy(),
                 8) == 8.52238967
    assert round(Water(3500, 700, massic_basis=True).entropy(),
                 7) == 10.1749996
    assert round(Water(30e6, 700, massic_basis=True).entropy(),
                 8) == 5.17540298

    #assert enthalpy
    assert round(Water(3500, 300, massic_basis=True).enthalpy(),
                 5) == 2549.91145
    assert round(Water(3500, 700, massic_basis=True).enthalpy(),
                 5) == 3335.68375
    assert round(Water(30e6, 700, massic_basis=True).enthalpy(),
                 5) == 2631.49474


    #assert cp
#    assert round(Water(3e6, 300).heat_capacity(),8) == 4.17301218
#    assert round(Water(80e6, 300).heat_capacity(),8) == 4.01008987
#    assert round(Water(3e6, 500).heat_capacity(),8) == 4.65580682

#    #assert enthalpy
#    assert round(Water(3e6, 300).enthalpy(),6) == 115.331273
#    assert round(Water(80e6, 300).enthalpy(),6) == 184.142828
#    assert round(Water(3e6, 500).enthalpy(),6) == 975.542239
#
#    #assert enthalpy
#    assert round(Water(3e6, 300).enthalpy(),6) == 115.331273
#    assert round(Water(80e6, 300).enthalpy(),6) == 184.142828
#    assert round(Water(3e6, 500).enthalpy(),6) == 975.542239

#region 3
    #assert specific volume
#    assert round(Water(3500, 300).specific_volume(),7) == 39.4913866
#    assert round(Water(3500, 700).specific_volume(),7) == 92.3015898
#    assert round(Water(30e6, 700).specific_volume(),11) == 0.00542946619
#
#    #assert internal energy
#    assert round(Water(3e6, 300).enthalpy(),6) == 115.331273
#    assert round(Water(80e6, 300).enthalpy(),6) == 184.142828
#    assert round(Water(3e6, 500).enthalpy(),6) == 975.542239
#
#    #assert enthropy
#    assert round(Water(3e6, 300).enthalpy(),6) == 115.331273
#    assert round(Water(80e6, 300).enthalpy(),6) == 184.142828
#    assert round(Water(3e6, 500).enthalpy(),6) == 975.542239

    #assert enthalpy
    assert round(Water(25.5837018e6, 650,
                       massic_basis=True).enthalpy(), 5) == 1863.43019
    assert round(Water(22.2930643e6, 650,
                       massic_basis=True).enthalpy(),
                 5) == round(2375.12401, 3)
    assert round(Water(78.3095639e6, 750,
                       massic_basis=True).enthalpy(), 5) == 2258.68845

    #assert cp
#    assert round(Water(3e6, 300).heat_capacity(),8) == 4.17301218
#    assert round(Water(80e6, 300).heat_capacity(),8) == 4.01008987
#    assert round(Water(3e6, 500).heat_capacity(),8) == 4.65580682

#    #assert enthalpy
#    assert round(Water(3e6, 300).enthalpy(),6) == 115.331273
#    assert round(Water(80e6, 300).enthalpy(),6) == 184.142828
#    assert round(Water(3e6, 500).enthalpy(),6) == 975.542239
#
#    #assert enthalpy
#    assert round(Water(3e6, 300).enthalpy(),6) == 115.331273
#    assert round(Water(80e6, 300).enthalpy(),6) == 184.142828
#    assert round(Water(3e6, 500).enthalpy(),6) == 975.542239

# region 4
    # There is no region 4; instead region 4 is the saturation line

# region 5
    #assert specific volume
#    assert round(Water(3500, 300).specific_volume(),7) == 39.4913866
#    assert round(Water(3500, 700).specific_volume(),7) == 92.3015898
#    assert round(Water(30e6, 700).specific_volume(),11) == 0.00542946619
#
#    #assert internal energy
#    assert round(Water(3e6, 300).enthalpy(),6) == 115.331273
#    assert round(Water(80e6, 300).enthalpy(),6) == 184.142828
#    assert round(Water(3e6, 500).enthalpy(),6) == 975.542239
#
#    #assert enthropy
#    assert round(Water(3e6, 300).enthalpy(),6) == 115.331273
#    assert round(Water(80e6, 300).enthalpy(),6) == 184.142828
#    assert round(Water(3e6, 500).enthalpy(),6) == 975.542239

    #assert enthalpy
    assert round(Water(0.5e6, 1500,
                       massic_basis=True).enthalpy(), 5) == 5219.76855
    assert round(Water(30e6, 1500,
                       massic_basis=True).enthalpy(), 5) == 5167.23514
    assert round(Water(30e6, 2000,
                       massic_basis=True).enthalpy(), 5) == 6571.22604

    #assert cp
#    assert round(Water(3e6, 300).heat_capacity(),8) == 4.17301218
#    assert round(Water(80e6, 300).heat_capacity(),8) == 4.01008987
#    assert round(Water(3e6, 500).heat_capacity(),8) == 4.65580682

#    #assert enthalpy
#    assert round(Water(3e6, 300).enthalpy(),6) == 115.331273
#    assert round(Water(80e6, 300).enthalpy(),6) == 184.142828
#    assert round(Water(3e6, 500).enthalpy(),6) == 975.542239
#
#    #assert enthalpy
#    assert round(Water(3e6, 300).enthalpy(),6) == 115.331273
#    assert round(Water(80e6, 300).enthalpy(),6) == 184.142828
#    assert round(Water(3e6, 500).enthalpy(),6) == 975.542239


# other tests
def triple_point_test():
    triple_temperature = 273.16
    triple_pressure = 611.657
    triple_water = Water(triple_pressure, triple_temperature)
    assert triple_water.internal_energy() < 1e-5
    assert triple_water.entropy() < 1e-5
