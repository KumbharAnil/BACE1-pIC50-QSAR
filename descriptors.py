
# Molecular descriptor calculations using PaDEL and RDKit

from padelpy import from_smiles
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors


REQUIRED_DESCRIPTORS = ["SIC1", "SpMin4_Bhm", "SC-5"]


LIPINSKI_RULES = {
    "MolWt": lambda value: value <= 500,
    "LogP": lambda value: value <= 5,
    "HBD": lambda value: value <= 5,
    "HBA": lambda value: value <= 10,
}


def extract_descriptors(desc_dict):
    """Extract the descriptors required by the QSAR equation."""

    missing = [
        descriptor
        for descriptor in REQUIRED_DESCRIPTORS
        if descriptor not in desc_dict
    ]

    if missing:
        raise KeyError(f"Missing descriptors: {missing}")

    return (
        float(desc_dict["SIC1"]),
        float(desc_dict["SpMin4_Bhm"]),
        float(desc_dict["SC-5"]),
    )


def calculate_padel_descriptors(smiles):
    """Calculate the three required descriptors for one SMILES."""

    descriptors = from_smiles(smiles, fingerprints=False)

    if isinstance(descriptors, list):
        descriptors = descriptors[0]

    return extract_descriptors(descriptors)


def compute_physicochemical(smiles):
    """Calculate molecular properties and Lipinski violations."""

    molecule = Chem.MolFromSmiles(smiles)

    if molecule is None:
        raise ValueError("Could not parse SMILES with RDKit.")

    mw = Descriptors.MolWt(molecule)
    logp = Crippen.MolLogP(molecule)
    tpsa = rdMolDescriptors.CalcTPSA(molecule)

    hbd = Lipinski.NumHDonors(molecule)
    hba = Lipinski.NumHAcceptors(molecule)

    rotatable_bonds = Lipinski.NumRotatableBonds(molecule)
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(molecule)
    molar_refractivity = Crippen.MolMR(molecule)

    values = {
        "MolWt": mw,
        "LogP": logp,
        "HBD": hbd,
        "HBA": hba,
    }

    violations = sum(
        not rule(values[prop])
        for prop, rule in LIPINSKI_RULES.items()
    )

    return {
        "MolWt": round(mw, 2),
        "LogP": round(logp, 3),
        "TPSA": round(tpsa, 2),
        "HBD": hbd,
        "HBA": hba,
        "RotatableBonds": rotatable_bonds,
        "AromaticRings": aromatic_rings,
        "MolarRefractivity": round(molar_refractivity, 3),
        "Lipinski_Violations": violations,
    }
