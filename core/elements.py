"""

chemkit.core.elements
======================

This module consists of functions on chemical elements in the periodic table.
"""

from chemkit.core.elementData import element_data


def get_molar_mass(element: str):
    """
    return molar mass of element
    """
    try:
        case_insensitive_element = element.title()
        return element_data[case_insensitive_element]["mass_number"]
    except ValueError:
        return f"Element not found: {element}"


def get_element_info(element: str):
    """
    return information on the element
    """
    try:
        case_insensitive_element = element.title()
        return element_data[case_insensitive_element]
    except ValueError:
        return f"Element not found: {element}"


__all__ = ["get_molar_mass", "get_element_info"]
