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


# ToDo: Move the functional group object to its own file for better readability
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
        "Carbonic Acid/Carbpnic Ester": "[CX3](=[OX1])(O)O",
        "Carbonic Acid/Carbonic Acid-Ester": "[CX3](=[OX1])([OX2])[OX2H,OX1H0-1]",
        "Carbonic Ester (carbonic acid diester)": "C[OX2][CX3](=[OX1])[OX2]C",
        "Carboxylic acid": "[CX3](=O)[OX2H1]",
        "Carboxylic acid/conjugate base": "[CX3](=O)[OX1H0-,OX2H1]",
        "Cyanamide": "[NX3][CX2]#[NX1]",
        "Ester": "[#6][CX3](=O)[OX2H0][#6]",
        "Ketone": "[#6][CX3](=O)[#6]",
        "Ether": "[OD2]([#6])[#6]",
        "Hydrogen": "[H]",
        "Enamine": "[NX3][CX3]=[CX3]",
        "Glycine": "N[CX4H2][CX3](=[OX1])[O,N]",
        "Proline": "N1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[O,N]",
        "Alanine side chain": "[CH3X4]",
        "Arginine side chain": (
            "[CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]"
        ),
        "Aspargine side chain": "[CH2X4][CX3](=[OX1])[NX3H2]",
        "Apartate": "[CH2X4][CX3](=[OX1])[OH0-,OH]",
        "Cysteine side chain": "[CH2X4][SX2H,SX1H0-]",
        "Glutamate": "[CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]",
        "Isoleucine side chain": "[CHX4]([CH3X4])[CH2X4][CH3X4]",
        "Leucine side chain": "[CH2X4][CHX4]([CH3X4])[CH3X4]",
        "Lysine side chain": "[CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]",
        "Thioamide": "[NX3][CX3]=[SX1]",
        "Valine side chain": "[CHX4]([CH3X4])[CH3X4]",
        "Alanine side chain": "[CH3X4]",
        "Azide group": "[$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]",
        "Azide ion": "[$([NX1-]=[NX2+]=[NX1-]),$([NX1]#[NX2+]-[NX1-2])]",
        "Hydrazine": "[NX3][NX3]",
        "Hydrazone": "[NX3][NX2]=[*]",
        "Iminium": "[NX3+]=[CX3]",
        "Unsubstituted dicarboximide": "[CX3](=[OX1])[NX3H][CX3](=[OX1])",
        "Substituted dicarboxmide": "[CX3](=[OX1])[NX3H0]([#6])[CX3](=[OX1])",
        "Dicarboxdiimide": (
            "[CX3](=[OX1])[NX3H0]([NX3H0]([CX3](=[OX1]))[CX3](=[OX1]))[CX3](=[OX1])"
        ),
        "Nitrate group": "[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]",
        "Nitrate Anion": "[$([OX1]=[NX3](=[OX1])[OX1-]),$([OX1]=[NX3+]([OX1-])[OX1-])]",
        "Nitrile": "[NX1]#[CX2]",
        "Isonitrile": "[CX1-]#[NX2+]",
        "Nitroso group": "[NX2]=[OX1]",
        "N-Oxide": "[$([#7+][OX1-]),$([#7v5]=[OX1]);!$([#7](~[O])~[O]);!$([#7]=[#7])]",
        "Hydroxyl in Carboxylic acid": "[OX2H][CX3]=[OX1]",
        "Enol": "[OX2H][#6X3]=[#6]",
        "Phenol": "[OX2H][cX3]:[c]",
        "Enol or Phenol": "[OX2H][$(C=C),$(cc)]",
        "Hydroxyl_acidic": "[$([OH]-*=[!#6])]",
        "Peroxide groups": "[OX2,OX1-][OX2,OX1-]",
        "Phosphoric acid groups": (
            "[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]"  # noqa E501
        ),
        "Phosphoric ester groups": (
            "[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]"  # noqa E501
        ),
        "Carbo-Thiocarboxylate": "[S-][CX3](=S)[#6]",
        "Carbo-Thioester": "S([#6])[CX3](=O)[#6]",
        "Thio/sulfide/disulfide sulfur": "[SX2]",
        "Thioamide": "[NX3][CX3]=[SX1]",
        "Sulfide": "[#16X2H0]",
        "Sulfinate": "[$([#16X3](=[OX1])[OX2H0]),$([#16X3+]([OX1-])[OX2H0])]",
        "Sulfinic acid": (
            "[$([#16X3](=[OX1])[OX2H,OX1H0-]),$([#16X3+]([OX1-])[OX2H,OX1H0-])]"
        ),
        "Carbo-azosulfone": "[SX4](C)(C)(=O)=N",
        "Sulfonamide": (
            "[$([SX4](=[OX1])(=[OX1])([!O])[NX3]),$([SX4+2]([OX1-])([OX1-])([!O])[NX3])]"  # noqa E501
        ),
        "Sulfenic acid": "[#16X2][OX2H,OX1H0-]",
        "Sulfenate": "[#16X2][OX2H0]",
        "Amino acid (generic)": "[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[O,N]",
        "Primary or secondary amine": "[NX3;H2,H1;!$(NC=O)]",
        "Primary amine": "[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6]",
        "Threonine side chain": "[CHX4]([CH3X4])[OX2H]",
        "Tryptophan side chain": (
            "[CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12"
        ),
        "Tyrosine side chain": (
            "[CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1"
        ),
        "Methionine side chain": "[CH2X4][CH2X4][SX2][CH3X4]",
        "Phenylalanine side chain": "[CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1",
        "Serine side chain": "[CH2X4][OX2H]",
    }

    mol = Chem.MolFromSmiles(smiles)

    detected = {}

    for name, smarts in functional_groups.items():
        pattern = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            detected[name] = matches
    return detected


__all__ = ["view_molecule", "get_molecular_formula", "get_molecular_similarity_factor"]
