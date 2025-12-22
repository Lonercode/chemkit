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

import py3Dmol
from rdkit import Chem
from rdkit.Chem import Draw


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
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            symbol = atom.GetSymbol()
            view.addLabel(
                symbol,
                {
                    "position": {
                        "x": mol.GetConformer().GetAtomPosition(idx).x,
                        "y": mol.GetConformer().GetAtomPosition(idx).y,
                        "z": mol.GetConformer().GetAtomPosition(idx).z,
                    },
                    "backgroundColor": "white",
                    "fontColor": "black",
                    "fontSize": 25,
                    "showBackground": False,
                    "fontOpacity": 1,
                },
            )
        view.setViewStyle(
            {
                "style": "outline",
                "color": "white",
                "width": 0.05,
            }
        )
        view.zoomTo()
        with open("molecule.html", "w") as f:
            f.write(view._make_html())
    else:
        print("Invalid molecule provided")
