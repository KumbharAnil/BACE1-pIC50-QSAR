# BACE1 QSAR model calculations

INACTIVE_CUTOFF = 4.0
MOST_ACTIVE_CUTOFF = 11.0


def compute_pic50(sic1, spmin4_bhm, sc5):
    """Calculate predicted pIC50 and IC50 in nanomolar units."""

    pIC50 = (
        -15.13463
        + (13.84656 * sic1)
        + (7.33758 * spmin4_bhm)
        + (1.90156 * sc5)
    )

    IC50_nM = 10 ** (9 - pIC50)

    return pIC50, IC50_nM


def classify_activity(pIC50):
    """Classify compounds using the existing activity thresholds."""

    if pIC50 < INACTIVE_CUTOFF:
        return "Inactive"

    if pIC50 > MOST_ACTIVE_CUTOFF:
        return "Most Active"

    return "Active"
