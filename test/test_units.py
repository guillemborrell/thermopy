# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 08:38:33 2015

@author: monteiro
"""
from thermopy.units import Pressure, Temperature


def test_units():
    """Test units."""
    Tkel = Temperature(273.15).unit('K')
    Tcel = Temperature(0).unit('C')
    assert Tkel == Tcel
    pressuremmhg = Pressure(760).unit('torr')
    pressureatm = Pressure(1).unit('atm')
    assert pressuremmhg == pressureatm
