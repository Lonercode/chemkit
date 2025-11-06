"""

chemkit.core.constants unit test
================================

Contains unit tests for chemkit.core.constants.py
"""

from chemkit.core import constants as c


def test_avogadro_constant():
    """Check that Avogadro's constant is accurate."""
    assert c.AVOGADRO == 6.02214076e23


def test_boltzmann_constant():
    """Check that Boltzmann's constant is accurate"""
    assert c.BOLTZMANN == 1.380651e-23


def test_faraday_constant():
    """Check that Faraday's constant is accurate"""
    assert c.FARADAY == 9.6485338e4


def test_gas_constant():
    """Check that Gas constant is accurate"""
    assert c.GAS_CONSTANT == 8.31447


def test_electron_mass_constant():
    """Check that electron mass is accurate"""
    assert c.ELECTRON_MASS == 9.109383e-28


def test_proton_mass_constant():
    """Check that proton mass is accurate"""
    assert c.PROTON_MASS == 1.6726217e-24


def test_neutron_mass_constant():
    """Check that Boltzmann's constant is accurate"""
    assert c.NEUTRON_MASS == 1.6749273e-24


def test_PI_constant():
    """Check that PI constant is accurate"""
    assert c.PI == 3.1415927


def test_planck_constant():
    """Check that Planck's constant is accurate"""
    assert c.PLANCK == 6.626069e-34


def test_light_speed_constant():
    """Check that spped of light constant is accurate"""
    assert c.SPEED_OF_LIGHT == 2.99792458e8


def test_compton_wavelength_electron_constant():
    """Check that the Compton wavelength of the electron is accurate"""
    assert c.COMPTON_WAVELENGTH_ELECTRON == 2.42631e-12


def test_bohr_radius_constant():
    """Check that Bohr radius is accurate"""
    assert c.BOHR_RADIUS == 5.29177e-11
