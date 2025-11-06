"""

chemkit.core.constants
======================

This module consists of chemical and physical constants
that would be used throughout the chemkit package.

Example
-------
>>> from chemkit.core.constants import AVOGADRO, GAS_CONSTANT
>>> GAS_CONSTANT / AVOGADRO
1.380649e-23
"""

from typing import Final

# Fundamental Physical Constants (SI Units)

AVOGADRO: Final[float] = 6.02214076e23  # mol^-1
BOLTZMANN: Final[float] = 1.380651e-23  # J.K^-1
FARADAY: Final[float] = 9.6485338e4  # C.mol^-1
GAS_CONSTANT: Final[float] = 8.31447  # J.molK^-1
ELECTRON_MASS: Final[float] = 9.109383e-28  # g
PROTON_MASS: Final[float] = 1.6726217e-24  # g
NEUTRON_MASS: Final[float] = 1.6749273e-24  # g
PI: Final[float] = 3.1415927
PLANCK: Final[float] = 6.626069e-34  # J.S
SPEED_OF_LIGHT: Final[float] = 2.99792458e8  # m.s^-1
COMPTON_WAVELENGTH_ELECTRON: Final[float] = 2.42631e-12  # m
BOHR_RADIUS: Final[float] = 5.29177e-11  # m

__all__ = [
    "AVOGADRO",
    "BOLTZMANN",
    "FARADAY",
    "GAS_CONSTANT",
    "ELECTRON_MASS",
    "PROTON_MASS",
    "NEUTRON_MASS",
    "PI",
    "SPEED_OF_LIGHT",
    "COMPTON_WAVELENGTH_ELECTRON",
    "BOHR_RADIUS",
]
