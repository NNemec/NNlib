from math import pi

# based upon:
# http://physics.nist.gov/cuu/Constants/index.html?/codata86.html

# numerical fundamental values
angstrom = 1.0
electron = 1.0         # electron-charge
eV = 1.0
hbar = 1.0
Kelvin = 1.0         # Kelvin

# derived values
nm = angstrom * 10
meter = angstrom * 1e10
Joule =  eV * 6.241506e+18
Coulomb = electron / 1.60217733e-19
bohr = 5.2918937910e-1 * angstrom

Volt = eV / electron

planck = 2*pi*hbar
second = planck / (6.6260755e-34 * Joule) # second

lightspeed = 2.99792458e+8 * meter / second # speed of light

kg = Joule * second**2 / meter**2
Ampere = Joule / (Volt * second)  # Ampere
Tesla = kg/(second**2 * Ampere)  # Tesla
kT = 1e3*Tesla

Phi_0 = planck / (2 * electron) # flux quantum
G_0 = 2*electron**2 / planck # conductance quantum

k_B = 8.617343e-5*eV/Kelvin   # Boltzmann

# basic quantities
a = 1.4226 * angstrom # basic length unit
#B_0 = 2 * Phi_0 / a^2
gamma = 2.66 * eV # basic coupling strength

m_e = 9.1093897e-31 * kg # electron mass
mu_B = electron*planck/(4*pi*m_e) # Bohr magneton
