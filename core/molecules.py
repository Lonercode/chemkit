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
from collections import defaultdict

import py3Dmol
from rdkit import Chem
from rdkit.Chem import DataStructs, Draw, rdFingerprintGenerator
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# ToDo: Work on symbols and metadata for molecule.


def view_molecule(dim: str, smiles: str = None, file: str = None):
    """View molecules"""
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
    """Obtain molecular formula for molecule based on smiles"""
    if file is None:
        mol = Chem.MolFromSmiles(smiles)
    else:
        mol = Chem.MolFromMolFile(file)

    formula = CalcMolFormula(mol)
    return formula


def get_molecular_similarity_factor(smiles1: str, smiles2: str):
    """Compare smiles for two molecules to determine degree of similarity"""
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    morgan_fngerprt = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2084)

    fngrprt1 = morgan_fngerprt.GetFingerprint(mol1)

    fngrprt2 = morgan_fngerprt.GetFingerprint(mol2)
    # check for similarity (the closer to 1, the greater the similarity)
    similarity = DataStructs.TanimotoSimilarity(fngrprt1, fngrprt2)
    return similarity


# ToDo: Add all functional groups avilable on
# https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
# in PR


def detect_functional_group(smiles: str):
    """Detect functional group(s) for molecules"""
    functional_groups = {
        "alkyl_carbon": "[CX4]",
        "allenic_carbon": "[$([CX2](=C)=C)]",
        "vinylic_carbon": "[$([CX3]=[CX3])]",
        "acetylenic_carbon": "[$([CX2]#C)]",
        "arene": "c",
        "carbony_group_low_specificity": "[CX3]=[OX1]",
        "carbonyl_group": "[$([CX3]=[OX1]),$([CX3+]-[OX1-])]",
        "carbonyl_with_carbon": "[CX3](=[OX1])C",
        "carbonyl_with_nitrogen": "[OX1]=CN",
        "carbonyl_with_oxygen": "[CX3](=[OX1])O",
        "acyl_halide": "[CX3](=[OX1])[F,Cl,Br,I]",
        "aldehyde": "[CX3H1](=O)[#6]",
        "anhydride": "[CX3](=[OX1])[OX2][CX3](=[OX1])",
        "amide": "[NX3][CX3](=[OX1])[#6]",
        "amidinium": "[NX3][CX3]=[NX3+]",
        "carbamate": "[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]",
        "carbamic_ester": "[NX3][CX3](=[OX1])[OX2H0]",
        "carbamic_acid": "[NX3,NX4+][CX3](=[OX1])[OX2H,OX1-]",
        "carboxylate_ion": "[CX3](=O)[O-]",
        "hydroxyl": "[OX2H]",
    }

    mol = Chem.MolFromSmiles(smiles)

    detected = {}

    for name, smarts in functional_groups.items():
        pattern = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            detected[name] = matches
    return detected


def element_count(smiles: str):
    """Get element count for each element in molecule"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")
    counts = defaultdict(int)

    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        counts[symbol] += 1
        counts["H"] += atom.GetTotalNumHs()

    return dict(counts)


def detect_degree_of_unsaturation(smiles: str):
    """Calculate degree of unsaturation"""
    counts = element_count(smiles)

    C = counts.get("C", 0)
    H = counts.get("H", 0)
    N = counts.get("N", 0)

    halogens = ["F", "Cl", "Br", "I"]

    X = sum(counts.get(x, 0) for x in halogens)

    du = (2 * C + 2 + N - H - X) / 2

    # A check since a negative degree of unsaturation would be nonsensical

    if du < 0:
        raise ValueError(" The degree of unsaturation is a negative value")

    return du


__all__ = ["view_molecule", "get_molecular_formula", "get_molecular_similarity_factor"]
