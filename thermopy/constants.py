u"""
Collection of physical constants and conversion factors.

Most constants are in SI units, so you can do
print('10 mile per minute is', 10*mile/minute, 'm/s or',
    10*mile/(minute*knot), 'knots')

The list is not meant to be comprehensive, but just a convenient list for
everyday use.
"""

from __future__ import absolute_import
import math as _math
# physical constants
# physical constants:
# (value, unit, precision)
cu_x_unit = (1.0020771e-13, 'm', 2.9e-20)
atomic_unit_of_time = (2.418884326505e-17, 's', 1.6e-28)
neutronproton_mass_ratio = (1.0013784187, '', 5.8e-10)
kilogramkelvin_relationship = (6.50965e+39, 'K', 1.1e+34)
atomic_mass_unitjoule_relationship = (1.4924179e-10, 'J', 2.6e-17)
hertzhartree_relationship = (1.519829846006e-16, 'E_h', 1e-27)
sackurtetrode_constant_1_k_101325_kpa = (-1.1648677, '', 4.4e-06)
atomic_unit_of_force = (8.2387225e-08, 'N', 1.4e-14)
tau_mass_energy_equivalent = (2.84705e-10, 'J', 4.6e-14)
muon_magnetic_moment = (-4.49044799e-26, 'J T^-1', 4e-33)
atomic_unit_of_1st_hyperpolarizablity = (3.20636151e-53, 'C^3 m^3 J^-2',
                                         2.8e-60)
electron_molar_mass = (5.4857990945e-07, 'kg mol^-1', 2.4e-16)
neutron_gyromagnetic_ratio_over_2_pi = (29.164695, 'MHz T^-1', 7.3e-06)
characteristic_impedance_of_vacuum = (376.73031346177066, 'ohm', 0.0)
atomic_unit_of_electric_field = (514220642000.0, 'V m^-1', 44000.0)
boltzmann_constant_in_inverse_meters_per_kelvin = (69.50356, 'm^-1 K^-1',
                                                   0.00012)
classical_electron_radius = (2.817940325e-15, 'm', 2.8e-23)
electronmuon_mass_ratio = (0.00483633167, '', 1.3e-10)
quantum_of_circulation_times_2 = (0.0007273895101, 'm^2 s^-1', 4.8e-12)
muonneutron_mass_ratio = (0.1124545175, '', 2.9e-09)
alpha_particle_mass_energy_equivalent_in_mev = (3727.37917, 'MeV', 0.00032)
molar_mass_of_carbon12 = (0.012, 'kg mol^-1', 0.0)
atomic_unit_of_electric_polarizablity = (1.648777274e-41, 'C^2 m^2 J^-1',
                                         1.6e-49)
atomic_mass_unithartree_relationship = (34231776.86, 'E_h', 0.23)
loschmidt_constant_27315_k_101325_kpa = (2.6867773e+25, 'm^-3', 4.7e+19)
molar_volume_of_ideal_gas_27315_k_100_kpa = (0.022710981, 'm^3 mol^-1', 4e-08)
proton_magnetic_moment_to_nuclear_magneton_ratio = (2.792847351, '', 2.8e-08)
atomic_unit_of_permittivity = (1.1126500560536183e-10, 'F m^-1', 0.0)
atomic_unit_of_2nd_hyperpolarizablity = (6.2353808e-65, 'C^4 m^4 J^-3',
                                         1.1e-71)
electron_voltinverse_meter_relationship = (806554.445, 'm^-1', 0.069)
inverse_meterjoule_relationship = (1.98644561e-25, 'J', 3.4e-32)
compton_wavelength_over_2_pi = (3.861592678e-13, 'm', 2.6e-21)
deuteronelectron_magnetic_moment_ratio = (-0.0004664345548, '', 5e-12)
neutronproton_magnetic_moment_ratio = (-0.68497934, '', 1.6e-07)
electron_voltatomic_mass_unit_relationship = (1.073544171e-09, 'u', 9.2e-17)
electron_mass_energy_equivalent = (8.1871047e-14, 'J', 1.4e-20)
electron_to_shielded_helion_magnetic_moment_ratio = (864.058255, '', 1e-05)
electronproton_magnetic_moment_ratio = (-658.2106862, '', 6.6e-06)
deuteronelectron_mass_ratio = (3670.4829652, '', 1.8e-06)
molar_mass_constant = (0.001, 'kg mol^-1', 0.0)
atomic_mass_unitkelvin_relationship = (10809527000000.0, 'K', 19000000.0)
lattice_parameter_of_silicon = (5.43102122e-10, 'm', 2e-17)
hertzkelvin_relationship = (4.7992374e-11, 'K', 8.4e-17)
natural_unit_of_action = (1.05457168e-34, 'J s', 1.8e-41)
tau_mass_energy_equivalent_in_mev = (1776.99, 'MeV', 0.29)
deuteron_magnetic_moment_to_bohr_magneton_ratio = (0.0004669754567, '', 5e-12)
standard_atmosphere = (101325.0, 'Pa', 0.0)
deuteronneutron_magnetic_moment_ratio = (-0.44820652, '', 1.1e-07)
wien_displacement_law_constant = (0.0028977685, 'm K', 5.1e-09)
electrondeuteron_mass_ratio = (0.00027244371095, '', 1.3e-13)
kelvinatomic_mass_unit_relationship = (9.251098e-14, 'u', 1.6e-19)
electron_magnetic_moment_to_nuclear_magneton_ratio = (-1838.28197107, '',
                                                      8.5e-07)
electron_volthertz_relationship = (241798940000000.0, 'Hz', 21000000.0)
helion_mass_energy_equivalent = (4.49953884e-10, 'J', 7.7e-17)
proton_mass_energy_equivalent = (1.50327743e-10, 'J', 2.6e-17)
hartreehertz_relationship = (6579683920721000.0, 'Hz', 44000.0)
muon_mass_in_u = (0.1134289264, 'u', 3e-09)
atomic_unit_of_electric_potential = (27.2113845, 'V', 2.3e-06)
atomic_mass_constant_energy_equivalent_in_mev = (931.494043, 'MeV', 8e-05)
two_hundred_twenty_lattice_spacing_of_silicon = (1.920155965e-10, 'm', 7e-18)
boltzmann_constant_in_evk = (8.617343e-05, 'eV K^-1', 1.5e-10)
electron_voltjoule_relationship = (1.60217653e-19, 'J', 1.4e-26)
atomic_unit_of_mass = (9.1093826e-31, 'kg', 1.6e-37)
conductance_quantum = (7.748091733e-05, 'S', 2.6e-13)
electron_volt = (1.60217653e-19, 'J', 1.4e-26)
neutron_g_factor = (-3.82608546, '', 9e-07)
electron_volthartree_relationship = (0.0367493245, 'E_h', 3.1e-09)
atomic_unit_of_current = (0.00662361782, 'A', 5.7e-10)
neutron_gyromagnetic_ratio = (183247183.0, 's^-1 T^-1', 46.0)
natural_unit_of_energy_in_mev = (0.510998918, 'MeV', 4.4e-08)
nuclear_magneton_in_kt = (0.00036582637, 'K T^-1', 6.4e-10)
kilogramatomic_mass_unit_relationship = (6.0221415e+26, 'u', 1e+20)
bohr_magneton_in_inverse_meters_per_tesla = (46.6864507, 'm^-1 T^-1', 4e-06)
natural_unit_of_momentum = (2.73092419e-22, 'kg m s^-1', 4.7e-29)
atomic_unit_of_action = (1.05457168e-34, 'J s', 1.8e-41)
electronneutron_mass_ratio = (0.00054386734481, '', 3.8e-13)
proton_mass_in_u = (1.00727646688, 'u', 1.3e-10)
rydberg_constant_times_c_in_hz = (3289841960360000.0, 'Hz', 22000.0)
hertzatomic_mass_unit_relationship = (4.439821667e-24, 'u', 3e-32)
neutron_to_shielded_proton_magnetic_moment_ratio = (-0.68499694, '', 1.6e-07)
atomic_unit_of_momentum = (1.99285166e-24, 'kg m s^-1', 3.4e-31)
hartreekilogram_relationship = (4.8508696e-35, 'kg', 8.3e-42)
kelvinjoule_relationship = (1.3806505e-23, 'J', 2.4e-29)
neutron_compton_wavelength_over_2_pi = (2.100194157e-16, 'm', 1.4e-24)
muon_g_factor = (-2.0023318396, '', 1.2e-09)
electronmuon_magnetic_moment_ratio = (206.7669894, '', 5.4e-06)
deuteron_rms_charge_radius = (2.1394e-15, 'm', 2.8e-18)
first_radiation_constant_for_spectral_radiance = (1.19104282e-16,
                                                  'W m^2 sr^-1', 2e-23)
atomic_unit_of_length = (5.291772108e-11, 'm', 1.8e-19)
tauelectron_mass_ratio = (3477.48, '', 0.57)
electron_mass = (9.1093826e-31, 'kg', 1.6e-37)
neutron_mass_energy_equivalent_in_mev = (939.56536, 'MeV', 8.1e-05)
shielded_helion_gyromagnetic_ratio_over_2_pi = (32.4341015, 'MHz T^-1',
                                                2.8e-06)
inverse_meterkilogram_relationship = (2.21021881e-42, 'kg', 3.8e-49)
alpha_particle_molar_mass = (0.004001506179149, 'kg mol^-1', 5.6e-14)
nuclear_magneton_in_mhzt = (7.62259371, 'MHz T^-1', 6.5e-07)
hartree_energy = (4.35974417e-18, 'J', 7.5e-25)
magnetic_constant = (1.2566370614359173e-06, 'N A^-2', 0.0)
electrontau_mass_ratio = (0.000287564, '', 4.7e-08)
proton_magnetic_moment = (1.41060671e-26, 'J T^-1', 1.2e-33)
deuteron_magnetic_moment = (4.33073482e-27, 'J T^-1', 3.8e-34)
rydberg_constant_times_hc_in_j = (2.17987209e-18, 'J', 3.7e-25)
tau_compton_wavelength_over_2_pi = (1.11046e-16, 'm', 1.8e-20)
quantum_of_circulation = (0.000363694755, 'm^2 s^-1', 2.4e-12)
planck_time = (5.39121e-44, 's', 4e-48)
thomson_cross_section = (6.65245873e-29, 'm^2', 1.3e-36)
second_radiation_constant = (0.014387752, 'm K', 2.5e-08)
neutron_mass_in_u = (1.0086649156, 'u', 5.5e-10)
tau_mass_in_u = (1.90768, 'u', 0.00031)
kilogramhartree_relationship = (2.06148605e+34, 'E_h', 3.5e+27)
proton_mass_energy_equivalent_in_mev = (938.272029, 'MeV', 8e-05)
protonneutron_mass_ratio = (0.99862347872, '', 5.8e-10)
atomic_mass_unitkilogram_relationship = (1.66053886e-27, 'kg', 2.8e-34)
hertzkilogram_relationship = (7.3724964e-51, 'kg', 1.3e-57)
neutron_mass = (1.67492728e-27, 'kg', 2.9e-34)
planck_constant = (6.6260693e-34, 'J s', 1.1e-40)
sackurtetrode_constant_1_k_100_kpa = (-1.1517047, '', 4.4e-06)
joulehartree_relationship = (2.29371257e+17, 'E_h', 39000000000.0)
electron_magnetic_moment = (-9.28476412e-24, 'J T^-1', 8e-31)
conventional_value_of_von_klitzing_constant = (25812.807, 'ohm', 0.0)
planck_constant_over_2_pi_in_ev_s = (6.58211915e-16, 'eV s', 5.6e-23)
electron_magnetic_moment_anomaly = (0.0011596521859, '', 3.8e-12)
proton_mass = (1.67262171e-27, 'kg', 2.9e-34)
natural_unit_of_time = (1.2880886677e-21, 's', 8.6e-30)
shielded_helion_gyromagnetic_ratio = (203789470.0, 's^-1 T^-1', 18.0)
hartreeelectron_volt_relationship = (27.2113845, 'eV', 2.3e-06)
neutron_magnetic_moment = (-9.6623645e-27, 'J T^-1', 2.4e-33)
shielded_helion_to_proton_magnetic_moment_ratio = (-0.761766562, '', 1.2e-08)
kilogramelectron_volt_relationship = (5.60958896e+35, 'eV', 4.8e+28)
deuteron_mass_in_u = (2.0135532127, 'u', 3.5e-10)
protontau_mass_ratio = (0.528012, '', 8.6e-05)
von_klitzing_constant = (25812.807449, 'ohm', 8.6e-05)
atomic_unit_of_charge_density = (1081202317000.0, 'C m^-3', 93000.0)
planck_constant_in_ev_s = (4.13566743e-15, 'eV s', 3.5e-22)
ideal_gas_constant = (8.314472, 'J mol^-1 K^-1', 1.5e-05)
ideal_gas_constant_massic_basis = (0.461526, 'J g^-1 K^-1', None)
finestructure_constant = (0.007297352568, '', 2.4e-11)
planck_constant_over_2_pi_times_c_in_mev_fm = (197.326968, 'MeV fm', 1.7e-05)
electron_g_factor = (-2.0023193043718, '', 7.5e-12)
tauproton_mass_ratio = (1.8939, '', 0.00031)
kelvinelectron_volt_relationship = (8.617343e-05, 'eV', 1.5e-10)
boltzmann_constant = (1.3806505e-23, 'J K^-1', 2.4e-29)
electron_voltkilogram_relationship = (1.78266181e-36, 'kg', 1.5e-43)
newtonian_constant_of_gravitation_over_hbar_c = (6.7087e-39, '(GeV/c^2)^-2',
                                                 1e-42)
boltzmann_constant_in_hzk = (20836644000.0, 'Hz K^-1', 36000.0)
atomic_unit_of_charge = (1.60217653e-19, 'C', 1.4e-26)
conventional_value_of_josephson_constant = (483597900000000.0, 'Hz V^-1', 0.0)
inverse_meterhartree_relationship = (4.55633525276e-08, 'E_h', 3e-19)
kilogramhertz_relationship = (1.35639266e+50, 'Hz', 2.3e+43)
proton_g_factor = (5.585694701, '', 5.6e-08)
tau_mass = (3.16777e-27, 'kg', 5.2e-31)
fermi_coupling_constant = (1.16639e-05, 'GeV^-2', 1e-10)
bohr_magneton = (9.27400949e-24, 'J T^-1', 8e-31)
deuteron_molar_mass = (0.0020135532127, 'kg mol^-1', 3.5e-13)
deuteron_mass_energy_equivalent = (3.00506285e-10, 'J', 5.1e-17)
kelvinkilogram_relationship = (1.5361808e-40, 'kg', 2.7e-46)
neutron_compton_wavelength = (1.3195909067e-15, 'm', 8.8e-24)
electronproton_mass_ratio = (0.00054461702173, '', 2.5e-13)
mo_x_unit = (1.00209966e-13, 'm', 5.3e-20)
helion_molar_mass = (0.0030149322434, 'kg mol^-1', 5.8e-12)
planck_constant_over_2_pi = (1.05457168e-34, 'J s', 1.8e-41)
kilograminverse_meter_relationship = (4.52443891e+41, 'm^-1', 7.7e+34)
rydberg_constant = (10973731.568525, 'm^-1', 7.3e-05)
nuclear_magneton_in_evt = (3.152451259e-08, 'eV T^-1', 2.1e-16)
bohr_magneton_in_evt = (5.788381804e-05, 'eV T^-1', 3.9e-13)
muonproton_mass_ratio = (0.1126095269, '', 2.9e-09)
standard_acceleration_of_gravity = (9.80665, 'm s^-2', 0.0)
electronneutron_magnetic_moment_ratio = (960.9205, '', 0.00023)
hartree_energy_in_ev = (27.2113845, 'eV', 2.3e-06)
hartreeinverse_meter_relationship = (21947463.13705, 'm^-1', 0.00015)
proton_gyromagnetic_ratio = (267522205.0, 's^-1 T^-1', 23.0)
shielded_helion_magnetic_moment = (-1.074553024e-26, 'J T^-1', 9.3e-34)
proton_rms_charge_radius = (8.75e-16, 'm', 6.8e-18)
natural_unit_of_action_in_ev_s = (6.58211915e-16, 'eV s', 5.6e-23)
deuteronproton_mass_ratio = (1.99900750082, '', 4.1e-10)
helion_mass_energy_equivalent_in_mev = (2808.39142, 'MeV', 0.00024)
nuclear_magneton_in_inverse_meters_per_tesla = (0.0254262358, 'm^-1 T^-1',
                                                2.2e-09)
atomic_unit_of_magnetic_dipole_moment = (1.8548019e-23, 'J T^-1', 1.6e-30)
inverse_of_conductance_quantum = (12906.403725, 'ohm', 4.3e-05)
electron_mass_in_u = (0.00054857990945, 'u', 2.4e-13)
kelvininverse_meter_relationship = (69.50356, 'm^-1', 0.00012)
neutron_magnetic_moment_to_bohr_magneton_ratio = (-0.00104187563, '', 2.5e-10)
joulekelvin_relationship = (7.242963e+22, 'K', 1.3e+17)
shielded_helion_magnetic_moment_to_bohr_magneton_ratio = (-0.001158671474, '',
                                                          1.4e-11)
natural_unit_of_momentum_in_mevc = (0.510998918, 'MeV/c', 4.4e-08)
magnetic_flux_quantum = (2.06783372e-15, 'Wb', 1.8e-22)
shielded_proton_magnetic_moment_to_nuclear_magneton_ratio = (2.792775604, '',
                                                             3e-08)
molar_volume_of_silicon = (1.20588382e-05, 'm^3 mol^-1', 2.4e-12)
tau_compton_wavelength = (6.9772e-16, 'm', 1.1e-19)
atomic_mass_unitinverse_meter_relationship = (751300660800000.0, 'm^-1',
                                              5000000.0)
muon_mass_energy_equivalent_in_mev = (105.6583692, 'MeV', 9.4e-06)
shielded_helion_magnetic_moment_to_nuclear_magneton_ratio = (-2.127497723, '',
                                                             2.5e-08)
muon_magnetic_moment_anomaly = (0.00116591981, '', 6.2e-10)
jouleinverse_meter_relationship = (5.0341172e+24, 'm^-1', 8.6e+17)
atomic_unit_of_electric_quadrupole_moment = (4.48655124e-40, 'C m^2', 3.9e-47)
kelvinhertz_relationship = (20836644000.0, 'Hz', 36000.0)
helion_mass_in_u = (3.0149322434, 'u', 5.8e-09)
inverse_meterkelvin_relationship = (0.014387752, 'K', 2.5e-08)
atomic_unit_of_electric_dipole_moment = (8.47835309e-30, 'C m', 7.3e-37)
avogadro_constant = (6.0221415e+23, 'mol^-1', 1e+17)
hartreekelvin_relationship = (315774.65, 'K', 0.55)
electric_constant = (8.854187817620389e-12, 'F m^-1', 0.0)
angstrom_star = (1.00001509e-10, 'm', 9e-17)
faraday_constant = (96485.3383, 'C mol^-1', 0.0083)
first_radiation_constant = (3.74177138e-16, 'W m^2', 6.4e-23)
electrondeuteron_magnetic_moment_ratio = (-2143.923493, '', 2.3e-05)
atomic_mass_unithertz_relationship = (2.252342718e+23, 'Hz',
                                      1500000000000000.0)
proton_charge_to_mass_quotient = (95788337.6, 'C kg^-1', 8.2)
elementary_charge_over_h = (241798940000000.0, 'A J^-1', 21000000.0)
muon_molar_mass = (0.0001134289264, 'kg mol^-1', 3e-12)
jouleelectron_volt_relationship = (6.24150947e+18, 'eV', 530000000000.0)
muonproton_magnetic_moment_ratio = (-3.183345118, '', 8.9e-08)
muon_compton_wavelength = (1.173444105e-14, 'm', 3e-22)
helion_mass = (5.00641214e-27, 'kg', 8.6e-34)
joulekilogram_relationship = (1.1126500560536185e-17, 'kg', 0.0)
electron_mass_energy_equivalent_in_mev = (0.510998918, 'MeV', 4.4e-08)
proton_magnetic_shielding_correction = (2.5689e-05, '', 1.5e-08)
electron_charge_to_mass_quotient = (-175882012000.0, 'C kg^-1', 15000.0)
atomic_unit_of_magnetic_flux_density = (235051.742, 'T', 0.02)
alpha_particle_mass_in_u = (4.001506179149, 'u', 5.6e-11)
atomic_unit_of_energy = (4.35974417e-18, 'J', 7.5e-25)
neutronelectron_mass_ratio = (1838.6836598, '', 1.3e-06)
taumuon_mass_ratio = (16.8183, '', 0.0027)
shielded_proton_magnetic_moment_to_bohr_magneton_ratio = (0.001520993132, '',
                                                          1.6e-11)
deuteron_mass_energy_equivalent_in_mev = (1875.61282, 'MeV', 0.00016)
atomic_mass_unitelectron_volt_relationship = (931494043.0, 'eV', 80.0)
natural_unit_of_length = (3.861592678e-13, 'm', 2.6e-21)
atomic_mass_constant = (1.66053886e-27, 'kg', 2.8e-34)
proton_molar_mass = (0.00100727646688, 'kg mol^-1', 1.3e-13)
electron_gyromagnetic_ratio = (176085974000.0, 's^-1 T^-1', 15000.0)
stefanboltzmann_constant = (5.6704e-08, 'W m^-2 K^-4', 4e-13)
tau_molar_mass = (0.00190768, 'kg mol^-1', 3.1e-07)
neutron_molar_mass = (0.0010086649156, 'kg mol^-1', 5.5e-13)
atomic_unit_of_velocity = (2187691.2633, 'm s^-1', 0.0073)
hartreeatomic_mass_unit_relationship = (2.921262323e-08, 'u', 1.9e-16)
electron_to_shielded_proton_magnetic_moment_ratio = (-658.2275956, '', 7.1e-06)
proton_compton_wavelength = (1.3214098555e-15, 'm', 8.8e-24)
shielded_proton_gyromagnetic_ratio = (267515333.0, 's^-1 T^-1', 23.0)
alpha_particle_mass_energy_equivalent = (5.9719194e-10, 'J', 1e-16)
electron_gyromagnetic_ratio_over_2_pi = (28024.9532, 'MHz T^-1', 0.0024)
inverse_finestructure_constant = (137.03599911, '', 4.6e-07)
protonneutron_magnetic_moment_ratio = (-1.45989805, '', 3.4e-07)
inverse_meteratomic_mass_unit_relationship = (1.3310250506e-15, 'u', 8.9e-24)
alpha_particleproton_mass_ratio = (3.97259968907, '', 5.2e-10)
proton_gyromagnetic_ratio_over_2_pi = (42.5774813, 'MHz T^-1', 3.7e-06)
elementary_charge = (1.60217653e-19, 'C', 1.4e-26)
electron_voltkelvin_relationship = (11604.505, 'K', 0.02)
jouleatomic_mass_unit_relationship = (6700536100.0, 'u', 1100.0)
unified_atomic_mass_unit = (1.66053886e-27, 'kg', 2.8e-34)
kelvinhartree_relationship = (3.1668153e-06, 'E_h', 5.5e-12)
neutron_mass_energy_equivalent = (1.50534957e-10, 'J', 2.6e-17)
bohr_magneton_in_hzt = (13996245800.0, 'Hz T^-1', 1200.0)
hertzjoule_relationship = (6.6260693e-34, 'J', 1.1e-40)
helionelectron_mass_ratio = (5495.885269, '', 1.1e-05)
alpha_particle_mass = (6.6446565e-27, 'kg', 1.1e-33)
newtonian_constant_of_gravitation = (6.6742e-11, 'm^3 kg^-1 s^-2', 1e-14)
atomic_unit_of_electric_field_gradient = (9.71736182e+21, 'V m^-2',
                                          830000000000000.0)
planck_mass = (2.17645e-08, 'kg', 1.6e-12)
electron_to_alpha_particle_mass_ratio = (0.000137093355575, '', 6.1e-14)
protonelectron_mass_ratio = (1836.15267261, '', 8.5e-07)
protonmuon_mass_ratio = (8.88024333, '', 2.3e-07)
inverse_meterhertz_relationship = (299792458.0, 'Hz', 0.0)
molar_planck_constant_times_c = (0.11962656572, 'J m mol^-1', 8e-10)
muonelectron_mass_ratio = (206.7682838, '', 5.4e-06)
inverse_meterelectron_volt_relationship = (1.23984191e-06, 'eV', 1.1e-13)
molar_planck_constant = (3.990312716e-10, 'J s mol^-1', 2.7e-18)
proton_magnetic_moment_to_bohr_magneton_ratio = (0.001521032206, '', 1.5e-11)
electron_magnetic_moment_to_bohr_magneton_ratio = (-1.0011596521859, '',
                                                   3.8e-12)
neutronelectron_magnetic_moment_ratio = (0.00104066882, '', 2.5e-10)
shielded_helion_to_shielded_proton_magnetic_moment_ratio = (-0.7617861313, '',
                                                            3.3e-09)
alpha_particleelectron_mass_ratio = (7294.2995363, '', 3.2e-06)
muon_mass_energy_equivalent = (1.6928336e-11, 'J', 2.9e-18)
muon_mass = (1.8835314e-28, 'kg', 3.3e-35)
shielded_proton_gyromagnetic_ratio_over_2_pi = (42.5763875, 'MHz T^-1',
                                                3.7e-06)
atomic_mass_constant_energy_equivalent = (1.4924179e-10, 'J', 2.6e-17)
speed_of_light_in_vacuum = (299792458.0, 'm s^-1', 0.0)
neutron_magnetic_moment_to_nuclear_magneton_ratio = (-1.91304273, '', 4.5e-07)
joulehertz_relationship = (1.50919037e+33, 'Hz', 2.6e+26)
atomic_unit_of_magnetizability = (7.8910366e-29, 'J T^-2', 1.3e-36)
tauneutron_mass_ratio = (1.89129, '', 0.00031)
deuteron_mass = (3.34358335e-27, 'kg', 5.7e-34)
proton_compton_wavelength_over_2_pi = (2.103089104e-16, 'm', 1.4e-24)
planck_length = (1.61624e-35, 'm', 1.2e-39)
nuclear_magneton = (5.05078343e-27, 'J T^-1', 4.3e-34)
muontau_mass_ratio = (0.0594592, '', 9.7e-06)
josephson_constant = (483597879000000.0, 'Hz V^-1', 41000000.0)
kilogramjoule_relationship = (8.987551787368176e+16, 'J', 0.0)
faraday_constant_for_conventional_electric_current = (96485.336, 'C_90 mol^-1',
                                                      0.016)
deuteronproton_magnetic_moment_ratio = (0.3070122084, '', 4.5e-09)
deuteron_magnetic_moment_to_nuclear_magneton_ratio = (0.8574382329, '',
                                                      9.2e-09)
bohr_magneton_in_kt = (0.6717131, 'K T^-1', 1.2e-06)
molar_volume_of_ideal_gas_27315_k_101325_kpa = (0.022413996, 'm^3 mol^-1',
                                                3.9e-08)
rydberg_constant_times_hc_in_ev = (13.6056923, 'eV', 1.2e-06)
natural_unit_of_velocity = (299792458.0, 'm s^-1', 0.0)
weak_mixing_angle = (0.22215, '', 0.00076)
hertzinverse_meter_relationship = (3.3356409519815204e-09, 'm^-1', 0.0)
hertzelectron_volt_relationship = (4.13566743e-15, 'eV', 3.5e-22)
compton_wavelength = (2.426310238e-12, 'm', 1.6e-20)
bohr_radius = (5.291772108e-11, 'm', 1.8e-19)
planck_temperature = (1.41679e+32, 'K', 1.1e+28)
hartreejoule_relationship = (4.35974417e-18, 'J', 7.5e-25)
neutrontau_mass_ratio = (0.52874, '', 8.6e-05)
muon_magnetic_moment_to_bohr_magneton_ratio = (-0.00484197045, '', 1.3e-10)
shielded_proton_magnetic_moment = (1.41057047e-26, 'J T^-1', 1.2e-33)
helionproton_mass_ratio = (2.9931526671, '', 5.8e-09)
muon_magnetic_moment_to_nuclear_magneton_ratio = (-8.89059698, '', 2.3e-07)
natural_unit_of_energy = (8.1871047e-14, 'J', 1.4e-20)
neutronmuon_mass_ratio = (8.89248402, '', 2.3e-07)
natural_unit_of_mass = (9.1093826e-31, 'kg', 1.6e-37)
muon_compton_wavelength_over_2_pi = (1.867594298e-15, 'm', 4.7e-23)


# mathematical constants
pi = _math.pi
golden = golden_ratio = (1 + _math.sqrt(5)) / 2

# SI prefixes
yotta = 1e24
zetta = 1e21
exa = 1e18
peta = 1e15
tera = 1e12
giga = 1e9
mega = 1e6
kilo = 1e3
hecto = 1e2
deka = 1e1
deci = 1e-1
centi = 1e-2
milli = 1e-3
micro = 1e-6
nano = 1e-9
pico = 1e-12
femto = 1e-15
atto = 1e-18
zepto = 1e-21

# binary prefixes
kibi = 2 ** 10
mebi = 2 ** 20
gibi = 2 ** 30
tebi = 2 ** 40
pebi = 2 ** 50
exbi = 2 ** 60
zebi = 2 ** 70
yobi = 2 ** 80

# weight in kg
gram = 1e-3
metric_ton = 1e3
grain = 64.79891e-6
lb = pound = 7000 * grain  # avoirdupois
oz = ounce = pound / 16
stone = 14 * pound
long_ton = 2240 * pound
short_ton = 2000 * pound

troy_ounce = 480 * grain  # only for metals / gems
troy_pound = 12 * troy_ounce
carat = 200e-6

# angle in rad
degree = pi / 180
arcmin = arcminute = degree / 60
arcsec = arcsecond = arcmin / 60

# time in second
minute = 60.0
hour = 60 * minute
day = 24 * hour
week = 7 * day
year = 365 * day
Julian_year = 365.25 * day

# length in meter
inch = 0.0254
foot = 12 * inch
yard = 3 * foot
mile = 1760 * yard
mil = inch / 1000
pt = point = inch / 72  # typography
survey_foot = 1200.0 / 3937
survey_mile = 5280 * survey_foot
nautical_mile = 1852.0
fermi = 1e-15
angstrom = 1e-10
micron = 1e-6
au = astronomical_unit = 149597870691.0
light_year = Julian_year * speed_of_light_in_vacuum[0]
parsec = au / arcsec

# pressure in pascal
atm = standard_atmosphere[0]
bar = 1e5
torr = mmHg = atm / 760
psi = pound * standard_acceleration_of_gravity[0] / (inch * inch)

# area in meter**2
hectare = 1e4
acre = 43560 * foot ** 2

# volume in meter**3
litre = liter = 1e-3
gallon = gallon_US = 231 * inch ** 3  # US
pint = gallon_US / 8
fluid_ounce = fluid_ounce_US = gallon_US / 128
bbl = barrel = 42 * gallon_US  # for oil

gallon_imp = 4.54609e-3  # uk
fluid_ounce_imp = gallon_imp / 160

# speed in meter per second
kmh = 1e3 / hour
mph = mile / hour
mach = speed_of_sound = 340.3
knot = nautical_mile / hour

# temperature in kelvin
zero_Celsius = 273.15
degree_Fahrenheit = 1 / 1.8  # only for differences

# energy in joule
eV = electron_volt = elementary_charge  # * 1 Volt
calorie = calorie_th = 4.184
calorie_IT = 4.1868
erg = 1e-7
Btu_th = pound * degree_Fahrenheit * calorie_th / gram
Btu = Btu_IT = pound * degree_Fahrenheit * calorie_IT / gram
ton_TNT = 1e9 * calorie_th
kWh = kilo * hour
# Wh = watt_hour

# power in watt
# hp = horsepower = 550 * foot * pound * g

# force in newton
dyn = dyne = 1e-5
lbf = pound_force = pound * standard_acceleration_of_gravity[0]
kgf = kilogram_force = standard_acceleration_of_gravity  # * 1 kg

# functions for conversions that are not linear


def C2K(C):
    """Convert Celcius to Kelvin"""
    return C + zero_Celsius


def K2C(K):
    """Convert Kelvin to Celcius"""
    return K - zero_Celsius


def F2C(F):
    """Convert Fahrenheit to Celcius"""
    return (F - 32) / 1.8


def C2F(C):
    """Convert Celcius to Fahrenheit"""
    return 1.8 * C + 32


def F2K(F):
    """Convert Fahrenheit to Kelvin"""
    return C2K(F2C(F))


def K2F(K):
    """Convert Kelvin to Fahrenheit"""
    return C2F(K2C(K))

# optics


def lambda2nu(lambda_):
    """Convert wavelength to optical frequency"""
    return speed_of_light_in_vacuum[0] / lambda_


def nu2lambda(nu):
    """Convert optical frequency to wavelength"""
    return speed_of_light_in_vacuum[0] / nu
