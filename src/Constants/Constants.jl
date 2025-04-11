module Constants

export CGRAV, CLIGHT, K_BOLTZ, AMU, CGAS, CRAD, SIGMA_SB, MSUN, RSUN, LSUN, MP, ME, SECYEAR, log10_e

# fundamental(ish) constants
"""
    CGRAV = 6.67430e-8

Gravitational constant in units of ``\\mathrm{g}^{-1}~\\mathrm{cm}^{3}~\\mathrm{s}^{-2}``

Source: [CODATA 2018](https://www.doi.org/10.1103/RevModPhys.93.025010)

Has a relative uncertainty of ``2.2\\times10^{-5}``
"""
const CGRAV = 6.67430e-8
    
"""
    CLIGHT = 2.99792458e10

Speed of light in units of ``\\mathrm{cm}~\\mathrm{s}^{-1}``.

Source [CODATA 2018](https://www.doi.org/10.1103/RevModPhys.93.025010)

Exact constant
"""
const CLIGHT = 2.99792458e10

"""
    PLANCK_H = 6.62607015e-27

Planck constant in units of ``\\mathrm{erg}~\\mathrm{s}``

Source: [CODATA 2018](https://www.doi.org/10.1103/RevModPhys.93.025010)

Exact constant
"""
const PLANCK_H = 6.62607015e-27          

"""
    HBAR = PLANCK_H / (2 * pi) = 1.0545718176461565e-27

Reduced Planck constant in units of ``\\mathrm{erg}~\\mathrm{s}``

Source: [CODATA 2018](https://www.doi.org/10.1103/RevModPhys.93.025010)

Exact constant
"""
const HBAR = PLANCK_H / (2 * pi)

"""
    QE = (CLIGHT/10)*1.602176634e-19 = 4.803204712570263e-10

Electron charge in esu==(g cm^3 s^-2)^(1/2) converted from the electron charge in Coulomb.

Source [CODATA 2018](https://www.doi.org/10.1103/RevModPhys.93.025010)

Exact constant
"""
const QE = (CLIGHT/10)*1.602176634e-19

"""
    AVO = 6.02214076e23

Avogadro's constant in units of ``\\mathrm{mol}^{-1}``.

Source [CODATA 2018](https://www.doi.org/10.1103/RevModPhys.93.025010)

Exact constant
"""
const AVO = 6.02214076e23

"""
    K_BOLTZ = 1.380649e-16

Boltzmann constant in units of ``\\mathrm{erg}^{-1}~\\mathrm{K}^{-1}``.

Source [CODATA 2018](https://www.doi.org/10.1103/RevModPhys.93.025010)

Exact constant
"""
const K_BOLTZ = 1.380649e-16

"""
    CGAS = K_BOLTZ * AVO = 8.31446261815324e7

Ideal gas constant in units of ``\\mathrm{erg}~\\mathrm{K}^{-1}~\\mathrm{mol}^{-1}``

Exact constant
"""
const CGAS = K_BOLTZ * AVO

"""
    AMU = 1.6605390666e-24

Atomic mass unit in units of ``\\mathrm{g}``

Source: [CODATA 2018](https://www.doi.org/10.1103/RevModPhys.93.025010)

Has a relative uncertainty of ``3\\times10^{-10}``
"""
const AMU = 1.6605390666e-24

"""
    MN = 1.6749274980e-24

Neutron mass in units of ``\\mathrm{g}``

Source: [CODATA 2018](https://www.doi.org/10.1103/RevModPhys.93.025010)

Has a relative uncertainty of ``5.7\\times10^{-10}``
"""
const MN = 1.6749274980e-24 

"""
    MP = 1.67262192369e-24

Neutron mass in units of ``\\mathrm{g}``

Source: [CODATA 2018](https://www.doi.org/10.1103/RevModPhys.93.025010)

Has a relative uncertainty of ``3.1\\times10^{-10}``
"""
const MP = 1.67262192369e-24

"""
    ME = 9.1093837015e-28

Electron mass in units of ``\\mathrm{g}``

Source: [CODATA 2018](https://www.doi.org/10.1103/RevModPhys.93.025010)

Has a relative uncertainty of ``3.0\\times10^{-10}``
"""
const ME = 9.1093837015e-28

"""
    RBOHR = HBAR^2 / (ME * QE^2) = 5.291772111941798e-9

Bohr radius in units of ``\\mathrm{cm}``
"""
const RBOHR = HBAR^2 / (ME * QE^2)         # Bohr radius (cm)

"""
    FINE = QE^2 / (HBAR * CLIGHT) = 0.007297352565305213

Fine structure constant
"""
const FINE = QE^2 / (HBAR * CLIGHT)          # fine structure constant

"""
    EV_TO_ERG = 1.602176634e-12

Energy of ``1~\\mathrm{eV}`` in  ``\\mathrm{erg}``

Source: [CODATA 2018](https://www.doi.org/10.1103/RevModPhys.93.025010)
"""
const EV_TO_ERG = 1.602176634e-12

"""
    MEV_TO_ERG = 1e6 * EV_TO_ERG = 1.602176634e-6

Energy of ``1~\\mathrm{MeV}`` in  ``\\mathrm{erg}``
"""
const MEV_TO_ERGS = 1e6 * EV_TO_ERG

"""
    SIGMA_SB = (2*π^5*K_BOLTZ^4)/(15*CLIGHT^2*PLANCK_H^3) = 5.67037441918443e-5

Stefan-Boltzmann constant in units of ``\\mathrm{erg}~\\mathrm{cm}^{-2}\\mathrm{K}^{-4}\\mathrm{s}^{-1}``
"""
const SIGMA_SB = (2*π^5*K_BOLTZ^4)/(15*CLIGHT^2*PLANCK_H^3)

"""
    CRAD = SIGMA_SB * 4 / CLIGHT = 7.565733250280006e-15

Radiation density constant in units of ``\\mathrm{erg}~\\mathrm{cm}^{-3}\\mathrm{K}^{-4}``
"""
const CRAD = SIGMA_SB * 4 / CLIGHT

"""
    G_TIMES_MSUN = 1.3271244e26

Product of the gravitational constant times the mass of the Sun. In units of ``\\mathrm{cm}^3~\\mathrm{s}^{-2}``.

Taken from the ["IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties"](https://doi.org/10.48550/arXiv.1510.07674).
"""
const G_TIMES_MSUN = 1.3271244e26

"""
    G_TIMES_MEARTH = 3.986004e20

Product of the gravitational constant times the mass of the Earth. In units of ``\\mathrm{cm}^3~\\mathrm{s}^{-2}``.

Taken from the ["IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties"](https://doi.org/10.48550/arXiv.1510.07674).
"""
const G_TIMES_MEARTH = 3.986004e20

"""
    G_TIMES_MJUPITER = 1.2668653e23

Product of the gravitational constant times the mass of Jupiter. In units of ``\\mathrm{cm}^3~\\mathrm{s}^{-2}``.

Taken from the ["IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties"](https://doi.org/10.48550/arXiv.1510.07674).
"""
const G_TIMES_MJUPITER = 1.2668653e23

"""
    MSUN = G_TIMES_MSUN/CGRAV = 1.9884098706980504e33

Mass of the Sun in ``\\mathrm{g}`` computed from the product GM.
"""
const MSUN = G_TIMES_MSUN/CGRAV

"""
    MEARTH = G_TIMES_MEARTH/CGRAV = 5.972167867791379e27

Mass of the Earth in ``\\mathrm{g}`` computed from the product GM.
"""
const MEARTH = G_TIMES_MEARTH/CGRAV

"""
    MJUPITER = G_TIMES_MJUPITER/CGRAV = 1.89812459733605e30

Mass of Jupiter in ``\\mathrm{g}`` computed from the product GM.
"""
const MJUPITER = G_TIMES_MJUPITER/CGRAV

"""
    RSUN = 6.957e10

Solar radius in ``\\mathrm{cm}``.

Taken from the ["IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties"](https://doi.org/10.48550/arXiv.1510.07674).
"""
const RSUN = 6.957e10

"""
    LSUN = 3.828e33 

Solar luminosity in ``\\mathrm{erg}~\\mathrm{s}^{-1}``.

Taken from the ["IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties"](https://doi.org/10.48550/arXiv.1510.07674).
"""
const LSUN = 3.828e33

"""
    DAYYEAR = 365.25e0

Number of days in a year
"""
const DAYYEAR = 365.25e0

"""
    SECYEAR = DAYYEAR * 24 * 3600 = 3.15576e7

Number of seconds in a year.
"""
const SECYEAR = DAYYEAR * 24 * 3600          # seconds per year

"""
    LY = CLIGHT*SECYEAR = 9.4607304725808e17

Lightyear in units of ``\\mathrm{cm}``
"""
const LY = CLIGHT*SECYEAR
const PC = 3.261633e0 * LY                 # parsec (cm)

"""
    AU = 1.49597870700e13

Astronomical unit in ``\\mathrm{cm}``.

From ["The IAU 2009 system of astronomical constants: the report of the IAU working group on numerical standards for Fundamental Astronomy"](https://doi.org/10.1007/s10569-011-9352-4).
"""
const AU = 1.49597870700e13


end
