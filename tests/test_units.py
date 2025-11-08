"""

chemkit.core.units unit test
================================

Contains unit tests for chemkit.core.units.py
"""

from chemkit.core import units as u

value = 25


def test_to_meter():
    """Check that conversion to meter works"""
    assert round((value / 100), 2) == u.to_meter(value, "cm")
    assert round((value / 1e-12), 2) == u.to_meter(value, "pm")


def test_to_joules():
    """Check that conversion to joules works"""
    assert round((value * 4.184), 2) == u.to_joules(value, "cal")
    assert round((value * 1.602e-19), 2) == u.to_joules(value, "ev")


def test_to_liter():
    """Check that conversion to joules works"""
    assert round((value * 4.184), 2) == u.to_joules(value, "cal")
    assert round((value * 1.602e-19), 2) == u.to_joules(value, "ev")
