module Constants

# Most of these are taken from MESA, would be good to update them as I took all
# this numbers quite some time ago

export CGRAV, CLIGHT, KERG, AMU, CGAS, CRAD, BOLTZ_SIGMA, MSUN, RSUN, LSUN, MP, ME, SECYEAR, log10_e

# fundamental(ish) constants
const CGRAV = 6.67428e-8                 # gravitational constant (g^-1 cm^3 s^-2)
const PLANCK_H = 6.62606896e-27          # Planck's constant (erg s)
const HBAR = PLANCK_H / (2 * pi)         # reduced Plank's constant (erg s)
const QE = 4.80320440e-10                # electron charge (esu == (g cm^3 s^-2)^(1/2))
const AVO = 6.02214179e23                # Avogadro's constant (mole^-1)
const CLIGHT = 2.99792458e10             # speed of light in vacuum (cm s^-1)
const KERG = 1.3806504e-16               # Boltzmann's constant (erg K^-1)
const CGAS = KERG * AVO                  # ideal gas constant; (erg K^-1)

# nuclear constants
const AMU = 1.660538782e-24              # atomic mass unit (g)
const MN = 1.6749286e-24                 # neutron mass (g)
const MP = 1.6726231e-24                 # proton mass (g)
const ME = 9.1093826e-28                 # (was 9.1093897d-28) electron mass (g)
const RBOHR = HBAR^2 / (ME * QE^2)       # Bohr radius (cm)
const FINE = QE^2 / (HBAR * CLIGHT)      # fine structure constant ≈ 1/137 (dimless)
const HION = 13.605698140e0              # hydrogen ionization energy (eV)
const EV_TO_ERG = 1.602176487e-12        # electron volts in one erg
const MEV_TO_ERGS = 1e6 * EV_TO_ERG      # MeVs in one erg
const MEV_AMU = MEV_TO_ERGS / AMU
const QCONV = MEV_TO_ERGS * AVO

# radiation constants
const BOLTZ_SIGMA = 5.670400e-5          # boltzmann's sigma = a*c/4 (erg cm^-2 K^-4 s^-1)
const CRAD = BOLTZ_SIGMA * 4 / CLIGHT    # = radiation density constant, a (erg cm^-3 K^-4)
                                         # Prad = crad * T^4 / 3, approx 7.5657e-15
const WIENLAM = PLANCK_H * CLIGHT / (KERG * 4.965114232e0)
const WIENFRE = 2.821439372E0 * KERG / PLANCK_H
const RHONUC = 2.342e14                  # density of nucleus (g cm^3)

# astronomical constants
# solar age, L, and R values from Bahcall et al, ApJ 618 (2005) 1049-1056.
const MSUN = 1.9892e33                   # solar mass (g) <<< gravitational mass, not baryonic
const RSUN = 6.9598e10                   # solar radius (cm)
const LSUN = 3.8418e33                   # solar luminosity (erg s^-1)
const AGESUN = 4.57e9                    # solar age (years)
const LY = 9.460528e17                   # light year (cm)
const PC = 3.261633e0 * LY               # parsec (cm)
const DAYYEAR = 365.25e0                 # days per year
const SECYEAR = DAYYEAR * 24 * 3600      # seconds per year
const TEFFSOL = 5777.0e0                 # solar teff (K)
const LOGGSOL = 4.4378893534131256e0     # With mesa's default msol, rsol and standard_cgrav
const MBOLSUN = 4.746                    # Bolometric magnitude of the Sun
const M_EARTH = 5.9764e27                # earth mass (g) = 3.004424e-6 Msun
const R_EARTH = 6.37e8                   # earth radius (cm)
const AU = 1.495978921e13                # astronomical unit (cm)
const M_JUPITER = 1.8986e30              # jupiter mass (g) = 0.954454d-3 Msun
const R_JUPITER = 6.9911e9               # jupiter mean radius (cm)
const SEMIMAJOR_AXIS_JUPITER = 7.7857e13 # jupiter semimajor axis (cm)

# mathematical constants
const log10_e = log10(exp(1))


end
