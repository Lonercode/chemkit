"""

chemkit.core.constants
======================

This module consists of mathematical equations in Chemistry.
"""

from typing import Dict

from chemkit.core.constants import GAS_CONSTANT, PLANCK
from chemkit.core.elements import get_molar_mass


def get_moles(mass: float, element: str):
    """Get number of moles of substance"""

    return mass / get_molar_mass(element)


def enthalpy(internal_energy: float, pressure: float, volume: float):
    """Solve and return the value for enthalpy change"""

    enthalpy = internal_energy - (pressure * volume)
    return enthalpy


def gibbs_free_energy(enthalpy: float, entropy: float, temperature: float):
    """Solve and return Gibbs's free energy with
    enthalpy_change in units of kilojoule,
    entropy change in units of joule per Kelvin and temperature in Kelvin.
    """

    gibbs_energy = enthalpy - (temperature * entropy)
    return gibbs_energy


def ideal_gas(
    quantity: str,
    volume: float = None,
    num_of_moles: float = None,
    pressure: float = None,
    temperature: float = None,
):
    """
    Solve for the ideal gas law
    """
    if quantity == "pressure":
        return (num_of_moles * GAS_CONSTANT * temperature) / volume
    elif quantity == "moles":
        return (pressure * volume) / (GAS_CONSTANT * temperature)
    elif quantity == "volume":
        return (num_of_moles * GAS_CONSTANT * temperature) / pressure
    elif quantity == "temperature":
        (pressure * volume) / (GAS_CONSTANT * num_of_moles)
    else:
        raise ValueError(f"Unsupported: {quantity}")


def helmholtz_free_energy(internal_energy: float, temperature: float, entropy: float):
    """Solve and return Helmholtz's free energy"""

    helmholtz_energy = internal_energy - (temperature * entropy)
    return helmholtz_energy


def general_form_rate_law(
    rate_constant: float, concentrations: Dict[str, float], orders: Dict[str, float]
):
    """
    Solve the general form rate law equation
    """

    rate = rate_constant
    for species, conc in concentrations.items():
        order = orders.get(species)
        rate *= conc**order
    return rate


def photon_energy(frequency: float):
    """
    Solve for the energy of a photon at a particular frequency
    """

    energy_of_photon = PLANCK * frequency
    return energy_of_photon


__all__ = [
    "get_moles",
    "enthalpy",
    "gibbs_free_energy",
    "ideal_gas",
    "helmholtz_free_energy",
    "general_form_rate_law",
    "photon_energy",
]
