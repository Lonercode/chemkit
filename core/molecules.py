"""

chemkit.core.molecules
======================

This module consists of functions dealing with molecules.

Example
-------
>>> from chemkit.core.molecules import view_molecule
>>> smile = "CCO"
>>> image3d = view_molecule(smile, "3d")
"""

import webbrowser

import py3Dmol
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# ToDo: Work on symbols and metadata for molecule.


def view_molecule(dim: str, smiles: str = None, file: str = None):
    if file is None:
        mol = Chem.MolFromSmiles(smiles)
    else:
        mol = Chem.MolFromMolFile(file)

    if dim == "2d":
        Chem.SanitizeMol(mol)
        view = Draw.MolToImage(mol)
        view.show()
    elif dim == "3d":
        Chem.SanitizeMol(mol)
        view = py3Dmol.view(width=800, height=800)
        view.addModel(Chem.MolToMolBlock(mol), "mol")
        view.setStyle({"stick": {}})
        view.zoomTo()
        with open("molecule.html", "w") as f:
            f.write(view._make_html())
        webbrowser.open_new("molecule.html")
    else:
        raise ValueError("Invalid molecule provided")


def get_molecular_formula(smiles: str = None, file: str = None):
    if file is None:
        mol = Chem.MolFromSmiles(smiles)
    else:
        mol = Chem.MolFromMolFile(file)

    formula = CalcMolFormula(mol)
    return formula


__all__ = ["view_molecule", "get_molecular_formula"]
