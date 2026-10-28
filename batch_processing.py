
"""
BACE1 QSAR Batch Processing
---------------------------
Processes individual compounds and batches of SMILES.
"""

from padelpy import from_smiles

from model import compute_pic50, classify_activity
from descriptors import (
    extract_descriptors,
    compute_physicochemical,
)


# Maximum number of compounds processed in one batch
CHUNK_SIZE = 50

# Sanity-check threshold for QSAR descriptors
SANITY_BOUND = 50.0

# Physicochemical properties included in the output
PHYSCHEM_COLUMNS = [
    "MolWt",
    "LogP",
    "TPSA",
    "HBD",
    "HBA",
    "RotatableBonds",
    "AromaticRings",
    "MolarRefractivity",
    "Lipinski_Violations",
]


def blank_physchem():
    """Return empty physicochemical values."""

    return {
        column: ""
        for column in PHYSCHEM_COLUMNS
    }


def score_one(mol_id, smiles, descriptors):

    try:
        physchem = compute_physicochemical(smiles)
        physchem_note = ""

    except Exception as error:
        physchem = blank_physchem()
        physchem_note = f"physchem failed: {error}"

    try:
        sic1, spmin4_bhm, sc5 = extract_descriptors(
            descriptors
        )

        pIC50, IC50_nM = compute_pic50(
            sic1,
            spmin4_bhm,
            sc5,
        )

        suspect = max(
            abs(sic1),
            abs(spmin4_bhm),
            abs(sc5),
        ) > SANITY_BOUND

        if suspect:
            status = (
                "SUSPECT: descriptor value implausibly large"
            )
        else:
            status = "ok"

        if physchem_note:
            status = f"{status} | {physchem_note}"

        return {
            "ID": mol_id,
            "SMILES": smiles,
            "SIC1": sic1,
            "SpMin4_Bhm": spmin4_bhm,
            "SC-5": sc5,
            "pIC50": round(pIC50, 4),
            "IC50_nM": round(IC50_nM, 4),
            "Activity": (
                "N/A"
                if suspect
                else classify_activity(pIC50)
            ),
            **physchem,
            "status": status,
        }

    except Exception as error:
        status = f"FAILED: {error}"

        if physchem_note:
            status = f"{status} | {physchem_note}"

        return {
            "ID": mol_id,
            "SMILES": smiles,
            "SIC1": "",
            "SpMin4_Bhm": "",
            "SC-5": "",
            "pIC50": "",
            "IC50_nM": "",
            "Activity": "N/A",
            **physchem,
            "status": status,
        }