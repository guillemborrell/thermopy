u"""
Calculates water and steam properties.

Classes:

    Water: water object.

References: [1] The International Association for the Properties of Water and
Steam. Revised Release on the IAPWS Industrial Formulation 1997 for the
Thermodynamic Properties of Water and Steam. August 2007.
"""

from thermopy.units import Pressure, Temperature, Enthalpy
from numpy import array, sum, sqrt, log
from thermopy.constants import ideal_gas_constant,                            \
    ideal_gas_constant_massic_basis
import scipy.optimize


class Water(object):
    """
    Taken from The International Association for the Properties of Water and
    Steam. Lucerne, Switzerland. August 2007. Revised Release on the IAPWS
    Industrial Formulation 1997 for the Thermodynamic Properties of Water and
    Steam. Reference document: IAPWS-IF97.

    All units in SI and default
    is molar basis."""

    def __init__(self, p, T, massic_basis=False):
        u"""Initializes a Water object."""

        self.p = Pressure(p)
        self.T = Temperature(T)
        # check if water is specified by IAPWS-IF97 for these values
        if self.T < Temperature(273.15) or self.T > Temperature(2273.15):
            raise ValueError('Temperature ' + str(T) +
                             ' out of range (273.15 - 2273.15K)')
        if self.p > Pressure(100).unit('MPa') or self.p < 0:
            raise ValueError('Pressure ' + str(p) +
                             ' out of range (0 - 100 MPa)')
        if self.T > 1073.15 and self.p > Pressure(50).unit('MPa'):
            raise ValueError('p or T value out of range')

        # adjust R to mass or molar basis
        if massic_basis is True:
            self.R = ideal_gas_constant_massic_basis[0]  # kJ/(kg K);
        elif massic_basis is False:
            self.R = ideal_gas_constant[0]

    # constants
    # ideal gas constant was already instantiated
    Tc = 647.096    # Critical point temperature (K)
    pc = Pressure(22.064).unit('MPa')     # Critical point pressure (MPa)
    rhoc = 322        # Critical point density kg/m3
    Tt = 273.16     # Triple point temperature (K)
    pt = 611.657    # Triple point pressure (Pa)
    ht = 0.611783    # Enthalpy at triple point (J/kg)

    def temperature_saturation(self, p=None):
        """Yields Tsat given a pressure p."""
        # module deals with pressure in MPa
        if p is None:
            p = self.p.MPa
        else:
            p = Pressure(p).MPa

        if p < Pressure(611.213).unit('Pa').MPa or p > self.pc:
            raise ValueError('Pressure out of range.')

        # table 34
        ni = array([1167.0521452767,
                    -724213.16703206,
                    -17.073846940092,
                    12020.82470247,
                    -3232555.0322333,
                    14.91510861353,
                    -4823.2657361591,
                    405113.40542057,
                    -0.23855557567849,
                    650.17534844798], dtype='d')

        beta = p ** 0.25

        E = 1 * beta ** 2 + ni[2] * beta + ni[5]
        F = ni[0] * beta ** 2 + ni[3] * beta + ni[6]
        G = ni[1] * beta ** 2 + ni[4] * beta + ni[7]
        D = 2 * G / (-F - (F**2 - 4 * E * G)**0.5)

        return Temperature((ni[9] + D - ((ni[9] + D) ** 2 - 4 *
                                         (ni[8] + ni[9] * D)) ** 0.5) * 0.5)

    def pressure_saturation(self, T=None):
        """Yields Psat given a temperature T"""
        if T is None:
            T = self.T
        else:
            T = Temperature(T)

        if T < 273.15 or T > self.Tc:
            raise ValueError('Temperature out of range.')

        # table 34
        ni = array([1167.0521452767,
                    -724213.16703206,
                    -17.073846940092,
                    12020.82470247,
                    -3232555.0322333,
                    14.91510861353,
                    -4823.2657361591,
                    405113.40542057,
                    -0.23855557567849,
                    650.17534844798], dtype='d')

        v = T + ni[8] / (T - ni[9])
        A = 1 * v ** 2 + ni[0] * v + ni[1]
        B = ni[2] * v ** 2 + ni[3] * v + ni[4]
        C = ni[5] * v ** 2 + ni[6] * v + ni[7]
        return Pressure((2 * C / (-B + (B ** 2 - 4 * A * C) ** 0.5)) ** 4
                        ).unit('MPa')

    def _is_in_region(self):
        """Finds a region for the (T, p) point (see IAPWS-IF97 for details).
        The usefulness of these regions are to divide the physical properties
        of water into different sets of coefficients and equations."""

        # for the 2 - 3 boundary line
        ni = array([348.05185628969,
                    -1.1671859879975,
                    0.0010192970039326,
                    572.54459862746,
                    13.91883977887], dtype='d')

        theta = self.T
        pressure23 = Pressure(ni[0] + ni[1] * theta + ni[2] *
                              theta ** 2).unit('MPa')
        # exceptional cases
        if self.T == self.Tt and self.p == self.pt:
            return 1
        # regular cases
        if self.T >= Temperature(273.15) and self.T <= Temperature(623.15):
            if self.T < self.temperature_saturation(self.p):
                return 1
            else:
                return 2
        elif self.T > Temperature(623.15) and self.T <= Temperature(1073.15):
            if self.p > pressure23:
                return 3
            else:
                return 2
        elif self.T >= Temperature(1073.15) and self.T <= Temperature(2273.15):
            return 5
        else:
            raise Exception('Cannot assign region to the parameters p = ' +
                            str(self.p) + 'T = ' + str(self.T) + 'given.')

    def _basic_equation1(self, value='gamma'):
        """Returns basic equation 1 and its derivatives, ex: 'gamma',
        'gamma_tau', 'gamma_tau_pi', etc."""
        # region1
        Ii = array([0, 0, 0, 0, 0, 0, 0, 0, 1,
                    1, 1, 1, 1, 1, 2, 2, 2, 2,
                    2, 3, 3, 3, 4, 4, 4, 5, 8,
                    8, 21, 23, 29, 30, 31, 32], dtype='int')
        Ji = array([-2, -1, 0, 1, 2, 3, 4, 5,
                    -9, -7, -1, 0, 1, 3, -3,
                    0, 1, 3, 17, -4, 0, 6, -5,
                    -2, 10, -8, -11, -6, -29,
                    -31, -38, -39, -40, -41], dtype='int')
        ni = array([0.14632971213167,
                    -0.84548187169114,
                    -3.756360367204,
                    3.3855169168385,
                    -0.95791963387872,
                    0.15772038513228,
                    -0.016616417199501,
                    0.00081214629983568,
                    0.00028319080123804,
                    -0.00060706301565874,
                    -0.018990068218419,
                    -0.032529748770505,
                    -0.021841717175414,
                    -5.283835796993e-05,
                    -0.00047184321073267,
                    -0.00030001780793026,
                    4.7661393906987e-05,
                    -4.4141845330846e-06,
                    -7.2694996297594e-16,
                    -3.1679644845054e-05,
                    -2.8270797985312e-06,
                    -8.5205128120103e-10,
                    -2.2425281908e-06,
                    -6.5171222895601e-07,
                    -1.4341729937924e-13,
                    -4.0516996860117e-07,
                    -1.2734301741641e-09,
                    -1.7424871230634e-10,
                    -6.8762131295531e-19,
                    1.4478307828521e-20,
                    2.6335781662795e-23,
                    -1.1947622640071e-23,
                    1.8228094581404e-24,
                    -9.3537087292458e-26], dtype='d')
        pi = self.p / Pressure(16.53).unit('MPa')
        tau = Temperature(1386) / self.T

        if value == 'gamma':
            return sum(ni * ((7.1 - pi) ** Ii) *
                       ((tau - 1.222) ** Ji))
        elif value == 'gamma_tau':
            return sum(ni * ((7.1 - pi) ** Ii) * Ji *
                       ((tau - 1.222) ** (Ji-1)))
        elif value == 'gamma_tau_tau':
            return sum(ni * ((7.1 - pi) ** Ii) * Ji * (Ji - 1) *
                       ((tau - 1.222) ** (Ji - 2)))
        elif value == 'gamma_pi':
            return -1 * sum(ni * Ii * ((7.1 - pi) ** (Ii - 1)) *
                            (tau - 1.222) ** Ji)
        elif value == 'gamma_pi_pi':
            return sum(ni * Ii * (Ii - 1) * (7.1 - pi) ** (Ii - 2) *
                       (tau - 1.222) ** Ji)
        elif value == 'gamma_pi_tau' or value == 'gamma_tau_pi':
            return -1 * sum(ni * Ii * (7.1 - pi) ** (Ii - 1) *
                            Ji * (tau - 1.222) ** (Ji - 1))
        else:
            raise Exception('Function not assigned in _basic_equation1()')

    def _basic_equation2(self, value='gamma'):
        """Returns equation 1 and its derivatives, ex: 'gamma', 'gamma_tau',
        'gamma_tau_pi', etc."""
        # region2; ideal part
        J0 = array([0, 1, -5, -4, -3, -2, -1,
                    2, 3], dtype='d')
        n0 = array([-9.6927686500217,
                    10.086655968018,
                    -0.005608791128302,
                    0.071452738081455,
                    -0.40710498223928,
                    1.4240819171444,
                    -4.383951131945,
                    -0.28408632460772,
                    0.021268463753307], dtype='d')
        # region2; real part
        Ii = array([1, 1, 1, 1, 1, 2, 2, 2, 2,
                    2, 3, 3, 3, 3, 3, 4, 4, 4,
                    5, 6, 6, 6, 7, 7, 7, 8, 8,
                    9, 10, 10, 10, 16, 16, 18,
                    20, 20, 20, 21, 22, 23,
                    24, 24, 24], dtype='d')
        Ji = array([0, 1, 2, 3, 6, 1, 2, 4, 7,
                    36, 0, 1, 3, 6, 35, 1, 2,
                    3, 7, 3, 16, 35, 0, 11,
                    25, 8, 36, 13, 4, 10, 14,
                    29, 50, 57, 20, 35, 48,
                    21, 53, 39, 26, 40, 58], dtype='d')
        ni = array([-0.0017731742473213,
                    -0.017834862292358,
                    -0.045996013696365,
                    -0.057581259083432,
                    -0.05032527872793,
                    -3.3032641670203e-05,
                    -0.00018948987516315,
                    -0.0039392777243355,
                    -0.043797295650573,
                    -2.6674547914087e-05,
                    2.0481737692309e-08,
                    4.3870667284435e-07,
                    -3.227767723857e-05,
                    -0.0015033924542148,
                    -0.040668253562649,
                    -7.8847309559367e-10,
                    1.2790717852285e-08,
                    4.8225372718507e-07,
                    2.2922076337661e-06,
                    -1.6714766451061e-11,
                    -0.0021171472321355,
                    -23.895741934104,
                    -5.905956432427e-18,
                    -1.2621808899101e-06,
                    -0.038946842435739,
                    1.1256211360459e-11,
                    -8.2311340897998,
                    1.9809712802088e-08,
                    1.0406965210174e-19,
                    -1.0234747095929e-13,
                    -1.0018179379511e-09,
                    -8.0882908646985e-11,
                    0.10693031879409,
                    -0.33662250574171,
                    8.9185845355421e-25,
                    3.0629316876232e-13,
                    -4.2002467698208e-06,
                    -5.9056029685639e-26,
                    3.7826947613457e-06,
                    -1.2768608934681e-15,
                    7.3087610595061e-29,
                    5.5414715350778e-17,
                    -9.436970724121e-07], dtype='d')

        pi = self.p / Pressure(1).unit('MPa')
        tau = 540 / self.T
        # arrays
        if value == 'gamma':
            return (self._basic_equation2('gamma_0') +
                    self._basic_equation2('gamma_r'))
        elif value == 'gamma_0':
            return log(pi) + sum(n0 * (tau ** J0))
        elif value == 'gamma_r':
            return sum(ni * (pi ** Ii) * (tau - 0.5) ** Ji)
        elif value == 'gamma_0_pi':
            return 1/pi
        elif value == 'gamme_0_pi_pi':
            return -1/(pi ** 2)
        elif value == 'gamma_0_tau':
            return 0 + sum(n0 * J0 * (tau ** (J0 - 1)))
        elif value == 'gamma_0_tau_tau':
            return 0 + sum(n0 * J0 * (J0 - 1) * tau ** (J0 - 2))
        elif value == 'gamma_0_pi_tau' or value == 'gamma_0_tau_pi':
            return 0
        elif value == 'gamma_r_pi':
            return sum(ni * Ii * (pi ** (Ii - 1)) * ((tau - 0.5) ** Ji))
        elif value == 'gamma_r_tau':
            return sum(ni * (pi ** Ii) * Ji * (tau - 0.5) ** (Ji - 1))
        elif value == 'gamma_r_pi_pi':
            return sum(ni * Ii * (Ii - 1) * (pi ** (Ii - 2))
                       (tau - 0.5) ** Ji)
        elif value == 'gamma_r_tau_tau':
            return sum(ni * (pi ** Ii) * Ji * (Ji - 1) * (tau - 0.5) **
                       (Ji - 2))
        elif value == 'gamma_pi_tau' or value == 'gamma_tau_pi':
            pass

    def _basic_equation3(self, value='PHI'):
        """Returns equation 3 and its derivatives, ex: 'gamma',
        'gamma_tau', 'gamma_tau_pi', etc."""
        # do not forget that the first term is not used
        Ii = array([-1, 0, 0, 0, 0, 0, 0, 0,
                    1, 1, 1, 1, 2, 2, 2, 2, 2,
                    2, 3, 3, 3, 3, 3, 4, 4, 4,
                    4, 5, 5, 5, 6, 6, 6, 7, 8,
                    9, 9, 10, 10, 11], dtype='d')
        # do not forget that the first term is not used
        Ji = array([-1, 0, 1, 2, 7, 10, 12,
                    23, 2, 6, 15, 17, 0, 2, 6,
                    7, 22, 26, 0, 2, 4, 16,
                    26, 0, 2, 4, 26, 1, 3, 26,
                    0, 2, 26, 2, 26, 2, 26, 0,
                    1, 26], dtype='d')
        ni = array([1.0658070028513,
                    -15.732845290239,
                    20.944396974307,
                    -7.6867707878716,
                    2.6185947787954,
                    -2.808078114862,
                    1.2053369696517,
                    -0.0084566812812502,
                    -1.2654315477714,
                    -1.1524407806681,
                    0.88521043984318,
                    -0.64207765181607,
                    0.38493460186671,
                    -0.85214708824206,
                    4.8972281541877,
                    -3.0502617256965,
                    0.039420536879154,
                    0.12558408424308,
                    -0.2799932969871,
                    1.389979956946,
                    -2.018991502357,
                    -0.0082147637173963,
                    -0.47596035734923,
                    0.0439840744735,
                    -0.44476435428739,
                    0.90572070719733,
                    0.70522450087967,
                    0.10770512626332,
                    -0.32913623258954,
                    -0.50871062041158,
                    -0.022175400873096,
                    0.094260751665092,
                    0.16436278447961,
                    -0.013503372241348,
                    -0.014834345352472,
                    0.00057922953628084,
                    0.0032308904703711,
                    8.0964802996215e-05,
                    -0.00016557679795037,
                    -4.4923899061815e-05], dtype='d')
        # it looks like the first term is never used
        # except in the PHI equations
        Ii = Ii[1:]
        Ji = Ji[1:]
        ni = ni[1:]
        n1 = 1.0658070028513

        tau = self.Tc / self.T
        # from equations 28 and p given as input,
        # calculate rho to be used in PHI

        def obj(x):
            return (self.p - 1000 * x * self.R * self.T * (x / self.rhoc) *
                    (n1 / (x / self.rhoc) + sum(
                        ni * Ii * (x / self.rhoc) ** (Ii - 1) * tau ** Ji)))
        for i in range(1, 580, 1):
            a = i
            b = a + 1
            if obj(a) * obj(b) < 0:
                break
        rho = scipy.optimize.bisect(obj, a, b, xtol=1e-10)
        delta = rho / self.rhoc

        # returns PHI and found rho as tuple (PHI, rho)
        if value == 'PHI':
            return (n1 * log(delta) + sum(ni * delta ** Ii * tau ** Ji),
                    rho)
        elif value == 'PHI_delta':
            return (n1 / delta + sum(ni * Ii * delta ** (Ii - 1) * tau ** Ji),
                    rho)
        elif value == 'PHI_delta_delta':
            return (-n1 / (delta ** 2) + sum(ni * Ii * (Ii - 1) * delta **
                                             (Ii - 2) * tau ** Ji),
                    rho)
        elif value == 'PHI_tau':
            return (sum(ni * delta ** Ii * Ji * tau ** (Ji - 1)),
                    rho)
        elif value == 'PHI_tau_tau':
            return (sum(ni * delta ** Ii * Ji * (Ji - 1) * tau ** (Ji - 2)),
                    rho)
        elif value == 'PHI_delta_tau' or value == 'PHI_tau_delta':
            return (sum(ni * Ii * delta ** (Ii - 1) * Ji * tau ** (Ji - 1)),
                    rho)

    # there is no region 4; region 4 is the saturation line of water/vapor

    def _basic_equation5(self, value='gamma'):
        """Returns equation 1 and its derivatives, ex: 'gamma',
        'gamma_tau', 'gamma_tau_pi', etc."""
        J0 = array([0, 1, -3, -2, -1, 2], dtype='d')
        n0 = array([-13.179983674201,
                    6.8540841634434,
                    -0.024805148933466,
                    0.36901534980333,
                    -3.1161318213925,
                    -0.32961626538917], dtype='d')
        Ii = array([1, 1, 1, 2, 2, 3], dtype='d')
        Ji = array([1, 2, 3, 3, 9, 7], dtype='d')
        ni = array([0.0015736404855259,
                    0.00090153761673944,
                    -0.0050270077677648,
                    2.2440037409485e-06,
                    -4.1163275453471e-06,
                    3.7919454822955e-08], dtype='d')
        pi = self.p / Pressure(1).unit('MPa')
        tau = Temperature(1000) / self.T

        if value == 'gamma':
            return (self._basic_equation5('gamma_0')
                    + self._basic_equation5('gamma_r'))
        elif value == 'gamma_0':
            return log(pi) + sum(n0 * (tau ** J0))
        elif value == 'gamma_r':
            return sum(ni * (pi ** Ii) * (tau ** Ji))
        elif value == 'gamma_0_pi':
            return 1 / pi
        elif value == 'gamma_0_pi_pi':
            return -1 / (pi ** 2)
        elif value == 'gamma_0_tau':
            return sum(n0 * J0 * (tau ** (J0 - 1)))
        elif value == 'gamma_0_tau_tau':
            sum(n0 * J0 * (J0 - 1) * (tau ** (J0 - 2)))
        elif value == 'gamma_0_tau_pi' or value == 'gamma_0_pi_tau':
            return 0 + 0
        elif value == 'gamma_r_pi':
            return sum(ni * Ii * (pi ** (Ii - 1)) * (tau ** Ji))
        elif value == 'gamma_r_tau':
            return sum(ni * (pi ** Ii) * Ji * (tau ** (Ji - 1)))
        elif value == 'gamma_r_pi_pi':
            return sum(ni * Ii * (Ii - 1) * (pi ** (Ii - 2)) * (tau ** Ji))
        elif value == 'gamma_r_tau_tau':
            return sum(ni * (pi ** Ii) * Ji * (Ji - 1) * (tau ** (Ji - 2)))
        elif value == 'gamma_r_tau_pi' or 'gamma_r_pi_tau':
            return sum(ni * Ii * (pi ** (Ii - 1)) * Ji * (tau ** (Ji - 1)))

    def gibbs_energy(self, p=None, T=None):
        """Returns the Gibbs Energy given p and T in kJ/kg."""
        if p is None:
            p = self.p
        else:
            p = Pressure(p)
        if T is None:
            T = self.T
        else:
            T = Temperature(T)

        if self._is_in_region() == 1:
            return self._basic_equation1('gamma') * self.R * T
        elif self._is_in_region() == 2:
            return self._basic_equation2('gamma') * self.R * T
        elif self._is_in_region() == 3:
            pass
        elif self._is_in_region() == 4:
            pass
        elif self._is_in_region() == 5:
            pass

    def specific_volume(self, p=None, T=None):
        """Returns the specific volume in m³/kg"""
        if p is None:
            p = self.p
        else:
            p = Pressure(p)
        if T is None:
            T = self.T
        else:
            T = Temperature(T)

        if self._is_in_region() == 1:
            pi = self.p / Pressure(16.53).unit('MPa')
            return (self._basic_equation1('gamma_pi') * self.R * self.T * pi
                    / self.p * 1e3)  # because of kJ in R
        elif self._is_in_region() == 2:
            pi = self.p / Pressure(1).unit('MPa')
            return (pi * (self._basic_equation2('gamma_0_pi') +
                          self._basic_equation2('gamma_r_pi')) * self.R
                    * self.T / self.p * 1e3)
        elif self._is_in_region() == 3:
            return None
        elif self._is_in_region() == 4:
            pass
        elif self._is_in_region() == 5:
            pass

        return -1

    def internal_energy(self, p=None, T=None):
        """NOT IMPLEMENTED"""
        if p is None:
            p = self.p
        else:
            p = Pressure(p)
        if T is None:
            T = self.T
        else:
            T = Temperature(T)

        if self._is_in_region() == 1:
            pi = self.p / Pressure(16.53).unit('MPa')
            tau = Temperature(1386) / self.T
            return (tau * self._basic_equation1('gamma_tau') - pi *
                    self._basic_equation1('gamma_pi')) * self.R * self.T

        elif self._is_in_region() == 2:
            pi = self.p / Pressure(1).unit('MPa')
            tau = 540 / self.T
            return (tau * (self._basic_equation2('gamma_0_tau') +
                           self._basic_equation2('gamma_r_tau')) - pi *
                    (self._basic_equation2('gamma_0_pi') +
                     self._basic_equation2('gamma_r_pi'))) * self.R * self.T
        elif self._is_in_region() == 3:
            pass
        elif self._is_in_region() == 4:
            pass
        elif self._is_in_region() == 5:
            pass

        return -1

    def entropy(self, p=None, T=None):
        """NOT IMPLEMENTED"""
        if p is None:
            p = self.p
        else:
            p = Pressure(p)
        if T is None:
            T = self.T
        else:
            T = Temperature(T)

        if self._is_in_region() == 1:
            pi = self.p / Pressure(16.53).unit('MPa')
            tau = Temperature(1386) / self.T
            return (tau * self._basic_equation1('gamma_tau') -
                    self._basic_equation1('gamma')) * self.R
        elif self._is_in_region() == 2:
            pi = self.p / Pressure(1).unit('MPa')
            tau = 540 / self.T
            return ((tau * (self._basic_equation2('gamma_0_tau') +
                            self._basic_equation2('gamma_r_tau')) -
                     (self._basic_equation2('gamma_0') +
                      self._basic_equation2('gamma_r'))) * self.R)
        elif self._is_in_region() == 3:
            pass
        elif self._is_in_region() == 4:
            pass
        elif self._is_in_region() == 5:
            pass

        return -1

    def enthalpy(self, p=None, T=None):
        """Returns the enthalpy given p and T in kJ/kg."""
        if p is None:
            p = self.p
        else:
            p = Pressure(p)
        if T is None:
            T = self.T
        else:
            T = Temperature(T)

        if self._is_in_region() == 1:
            return Enthalpy(self._basic_equation1('gamma_tau')
                            * 1386 * self.R)

        elif self._is_in_region() == 2:
            return Enthalpy(540 * self.R * (
                            self._basic_equation2('gamma_0_tau') +
                            self._basic_equation2('gamma_r_tau')))

        elif self._is_in_region() == 3:
            # for region 3 rho needs to be solved numerically
            reg3_solved = self._basic_equation3('PHI_delta')
            delta = reg3_solved[1] / self.rhoc
            tau = self.Tc / self.T
            return Enthalpy((tau * self._basic_equation3('PHI_tau')[0]
                             + delta * reg3_solved[0]) * self.R * self.T)

        elif self._is_in_region() == 4:
            return None

        elif self._is_in_region() == 5:
            pi = self.p / Pressure(1).unit('MPa')
            tau = Temperature(1000) / self.T
            return (tau * (self._basic_equation5('gamma_0_tau') +
                           self._basic_equation5('gamma_r_tau')) * self.R *
                    self.T)

    def heat_capacity(self, p=None, T=None):
        """Returns the isobaric heat capacity given p and T in kJ/(kg K)."""
        if p is None:
            p = self.p
        else:
            p = Pressure(p)
        if T is None:
            T = self.T
        else:
            T = Temperature(T)

        if self._is_in_region() == 1:
            tau = 1386 / T
            return (-1 * tau * tau * self.R *
                    self._basic_equation1('gamma_tau_tau'))
        elif self._is_in_region() == 2:
            pass
        elif self._is_in_region() == 3:
            delta = rho / self.rhoc
            tau = self.Tc / self.T
            return (- tau ** 2 / self._basic_equation3('PHI_tau_tau') +
                    (delta * self._basic_equation3('PHI_delta') -
                     delta * tau * self._basic_equation3('PHI_tau_delta'))
                    ** 2 / (2 * delta * self._basic_equation3('PHI_delta') +
                            delta ** 2 * self._basic_equation3(
                                'PHI_delta_delta'))
                    ) * self.R
        elif self._is_in_region() == 4:
            pass
        elif self._is_in_region() == 5:
            pass

        return -1

    def heat_capacity_constant_v(self, p=None, T=None):
        """NOT IMPLEMENTED"""
        if p is None:
            p = self.p
        else:
            p = Pressure(p)
        if T is None:
            T = self.T
        else:
            T = Temperature(T)

        if self._is_in_region() == 1:
            pi = self.p / Pressure(16.53).unit('MPa')
            tau = Temperature(1386) / self.T
            return (- (tau ** 2) * self._basic_equation1('gamma_tau_tau') +
                    ((self._basic_equation1('gamma_pi') - tau *
                      self._basic_equation1('gamma_pi_tau')) ** 2) /
                    self._basic_equation1('gamma_pi_pi')) * self.R
        elif self._is_in_region() == 2:
            pass
        elif self._is_in_region() == 3:
            pass
        elif self._is_in_region() == 4:
            pass
        elif self._is_in_region() == 5:
            pass

        return -1

    def speed_of_sound(self, p=None, T=None):
        """NOT IMPLEMENTED"""
        if p is None:
            p = self.p
        else:
            p = Pressure(p)
        if T is None:
            T = self.T
        else:
            T = Temperature(T)

        if self._is_in_region() == 1:
            pi = self.p / Pressure(16.53).unit('MPa')
            tau = Temperature(1386) / self.T
            inside_frac = (((self._basic_equation1('gamma_pi') - tau *
                             self._basic_equation1('gamma_tau_pi')) ** 2) /
                           ((tau ** 2) *
                            self._basic_equation1('gamma_tau_tau')
                            ))
            return sqrt(((self._basic_equation1('gamma_pi') ** 2) / (
                          inside_frac - self._basic_equation1('gamma_pi_pi'))
                         ) * self.R * self.T * 1000)

        elif self._is_in_region() == 2:
            pass
        elif self._is_in_region() == 3:
            pass
        elif self._is_in_region() == 4:
            pass
        elif self._is_in_region() == 5:
            pass
        return -1
