from math import pi

# based upon:
# http://physics.nist.gov/cuu/Constants/index.html?/codata86.html

second = 1.0
meter = 1.0
kg = 1.0
Ampere = 1.0
Kelvin = 1.0

###################

Newton = meter * kg / second**2
Joule = Newton * meter
Watt = Joule / second
Coulomb = Ampere * second
Volt = Watt / Ampere
Ohm = Volt / Ampere
Tesla = kg/(second**2 * Ampere)

###################

nm = 1e-9 * meter
Angstrom = 1e-10 * meter

electron = 1.60217733e-19 * Coulomb
eV = electron * Volt

bohr = 5.2918937910e-2 * nm

planck = 6.6260755e-34 * Joule * second
hbar = planck/(2*pi)

lightspeed = 2.99792458e+8 * meter / second # speed of light

k_B = 8.617343e-5 * eV / Kelvin   # Boltzmann
Phi_0 = planck / electron # flux quantum
G_0 = 2 * electron**2 / planck # conductance quantum

m_e = 9.1093897e-31 * kg # electron mass

mu_B = electron*planck/(4*pi*m_e) # Bohr magneton

rydberg = 13.6056923 * eV

mol = 6.0221415e+23 # Avogadro constant
