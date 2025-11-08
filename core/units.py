"""

chemkit.core.units
==================

This module consists of units and corresponding conversion
functions that would be used throughout the chemkit package.

Example
-------
>>> from chemkit.core.units import to_meter
>>> length = 1200
>>> length = to_meter(length, "cm")
12
"""

from typing import Final

# Conversion factors
# length
CM_TO_METER: Final[float] = 100.0
PM_TO_METER: Final[float] = 1e-12
# volume
ML_TO_LITER: Final[float] = 1000.0
# pressure
ATM_TO_TORRES: Final[float] = 760.0
# energy
CAL_TO_JOULES: Final[float] = 4.184
EV_TO_JOULES: Final[float] = 1.602e-19

# Functions that implement conversions


def to_meter(value: float, unit: str):
    """Convert value from unit to meter"""
    if unit == "cm":
        converted = round((value / CM_TO_METER), 2)
    elif unit == "pm":
        converted = round((value / PM_TO_METER), 2)
    else:
        raise ValueError(f"Unsupported unit: {unit}")
    return converted


def to_joules(value: float, unit: str):
    """Convert value from unit to joules"""
    if unit == "cal":  # Convert from calories
        converted = round((value * CAL_TO_JOULES), 2)
    elif unit == "ev":  # Convert from electronvolts
        converted = round((value * EV_TO_JOULES), 2)
    else:
        raise ValueError(f"Unsupported unit: {unit}")
    return converted


def to_liter(value: float, unit: str):
    """Convert value from unit to liter"""
    if unit == "ml":  # Convert from calories
        converted = round((value / ML_TO_LITER), 2)
    else:
        raise ValueError(f"Unsupported unit: {unit}")
    return converted


__all__ = [
    "CM_TO_METER",
    "PM_TO_METER",
    "ML_TO_LITER",
    "ATM_TO_TORRES",
    "CAL_TO_JOULES",
    "EV_TO_JOULES",
    "to_meter",
    "to_joules",
    "to_liter",
]
