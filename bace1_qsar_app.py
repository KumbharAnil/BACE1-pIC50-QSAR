"""
BACE1 QSAR pIC50 Prediction & Activity Classification Tool
------------------------------------------------------------
Streamlit app. The QSAR equation, descriptor requirements, sanity-check
bound, and chunked/fallback PaDEL calling strategy are preserved exactly
from the original Flask script — only the interface layer (Streamlit UI,
charts, activity classification) is new.

Activity classification (based on predicted pIC50):
    pIC50 <  4                -> INACTIVE
    4  <= pIC50 <= 11          -> ACTIVE
    pIC50 > 11                 -> MOST ACTIVE

Also computes standard physicochemical / Lipinski "Rule of Five"
descriptors for each SMILES (via RDKit) alongside the QSAR prediction.

Run with:
    streamlit run bace1_qsar_app.py
"""

import io
import csv
import traceback

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from padelpy import from_smiles
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
# ----------------------------------------------------------------------
# Page config & theme
# ----------------------------------------------------------------------

st.set_page_config(
    page_title="BACE1 QSAR Predictor",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

REQUIRED_DESCRIPTORS = ["SIC1", "SpMin4_Bhm", "SC-5"]
CHUNK_SIZE = 50

# SIC1, SpMin4_Bhm and SC-5 are small topological/information-content indices —
# in practice they sit well under this bound. A value outside it almost always
# means the wrong column got fed into PaDEL (e.g. an ID instead of a SMILES),
# not a real descriptor. Flag it instead of silently reporting "ok".
SANITY_BOUND = 50.0

# Activity classification cutoffs (on predicted pIC50)
INACTIVE_CUTOFF = 4.0
MOST_ACTIVE_CUTOFF = 11.0

# Physicochemical descriptor columns, in display order
PHYSCHEM_COLUMNS = [
    "MolWt", "LogP", "TPSA", "HBD", "HBA",
    "RotatableBonds", "AromaticRings", "MolarRefractivity", "Lipinski_Violations",
]

# Lipinski "Rule of Five" thresholds — one violation per rule broken
LIPINSKI_RULES = {
    "MolWt": (lambda v: v <= 500, "MW ≤ 500 Da"),
    "LogP": (lambda v: v <= 5, "LogP ≤ 5"),
    "HBD": (lambda v: v <= 5, "H-bond donors ≤ 5"),
    "HBA": (lambda v: v <= 10, "H-bond acceptors ≤ 10"),
}

ACTIVITY_COLORS = {
    "Inactive": "#d62728",      # red
    "Active": "#f2c744",        # amber
    "Most Active": "#2ca02c",   # green
    "N/A": "#7f7f7f",           # grey, used for failed/suspect rows
}
ACTIVITY_ORDER = ["Inactive", "Active", "Most Active", "N/A"]


# ----------------------------------------------------------------------
# Core scoring logic (unchanged from the original script)
# ----------------------------------------------------------------------

def compute_pic50(sic1, spmin4_bhm, sc5):
    pIC50 = (
        -15.13463
        + (13.84656 * sic1)
        + (7.33758 * spmin4_bhm)
        + (1.90156 * sc5)
    )
    IC50_nM = (10 ** (-pIC50)) * 1_000_000_000
    return pIC50, IC50_nM


def extract_descriptors(desc_dict):
    missing = [d for d in REQUIRED_DESCRIPTORS if d not in desc_dict]
    if missing:
        raise KeyError(f"missing descriptor(s): {missing}")
    return (
        float(desc_dict["SIC1"]),
        float(desc_dict["SpMin4_Bhm"]),
        float(desc_dict["SC-5"]),
    )


def compute_physicochemical(smi):
    """Standard physicochemical / Lipinski Ro5 descriptors for one SMILES,
    computed with RDKit. Returns a dict, or raises ValueError if the
    SMILES can't be parsed."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError("could not parse SMILES with RDKit")

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rot_bonds = Lipinski.NumRotatableBonds(mol)
    arom_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    mr = Crippen.MolMR(mol)

    values = {"MolWt": mw, "LogP": logp, "HBD": hbd, "HBA": hba}
    violations = sum(1 for prop, (rule, _) in LIPINSKI_RULES.items() if not rule(values[prop]))

    return {
        "MolWt": round(mw, 2),
        "LogP": round(logp, 3),
        "TPSA": round(tpsa, 2),
        "HBD": hbd,
        "HBA": hba,
        "RotatableBonds": rot_bonds,
        "AromaticRings": arom_rings,
        "MolarRefractivity": round(mr, 3),
        "Lipinski_Violations": violations,
    }


def classify_activity(pIC50):
    """Bucket a predicted pIC50 into the BACE1 activity classes."""
    if pIC50 < INACTIVE_CUTOFF:
        return "Inactive"
    if pIC50 > MOST_ACTIVE_CUTOFF:
        return "Most Active"
    return "Active"


def _blank_physchem():
    return {col: "" for col in PHYSCHEM_COLUMNS}


def score_one(mol_id, smi, desc):
    # Physicochemical properties are computed independently of the QSAR
    # descriptors/PaDEL — a SMILES that fails PaDEL can still get RDKit
    # properties, and vice versa, so failures here don't abort the row.
    try:
        physchem = compute_physicochemical(smi)
        physchem_note = ""
    except Exception as e:
        physchem = _blank_physchem()
        physchem_note = f"physchem failed: {e}"

    try:
        sic1, spmin4_bhm, sc5 = extract_descriptors(desc)
        pIC50, IC50_nM = compute_pic50(sic1, spmin4_bhm, sc5)
        suspect = max(abs(sic1), abs(spmin4_bhm), abs(sc5)) > SANITY_BOUND
        status = (
            "SUSPECT: descriptor value implausibly large — "
            "check ID/SMILES column mapping"
            if suspect else "ok"
        )
        if physchem_note:
            status = f"{status} | {physchem_note}"
        row = {
            "ID": mol_id,
            "SMILES": smi,
            "SIC1": sic1,
            "SpMin4_Bhm": spmin4_bhm,
            "SC-5": sc5,
            "pIC50": round(pIC50, 4),
            "IC50_nM": round(IC50_nM, 4),
            "Activity": "N/A" if suspect else classify_activity(pIC50),
            **physchem,
            "status": status,
        }
        return row
    except Exception as e:
        status = f"FAILED: {e}"
        if physchem_note:
            status = f"{status} | {physchem_note}"
        return {
            "ID": mol_id, "SMILES": smi, "SIC1": "", "SpMin4_Bhm": "",
            "SC-5": "", "pIC50": "", "IC50_nM": "", "Activity": "N/A",
            **physchem,
            "status": status,
        }


def predict_single(mol_id, smi):
    desc = from_smiles(smi, fingerprints=False)
    if isinstance(desc, list):
        desc = desc[0]
    return score_one(mol_id, smi, desc)


def process_chunk(chunk):
    """chunk = list of (mol_id, smi). Try a batched PaDEL call first,
    fall back to one-by-one so a single bad SMILES doesn't fail the batch."""
    results = []
    ids = [m for m, _ in chunk]
    smiles_list = [s for _, s in chunk]
    try:
        desc_list = from_smiles(smiles_list, fingerprints=False)
        if not isinstance(desc_list, list):
            desc_list = [desc_list]
        if len(desc_list) != len(chunk):
            raise ValueError("descriptor count mismatch — falling back")
        for mol_id, smi, desc in zip(ids, smiles_list, desc_list):
            results.append(score_one(mol_id, smi, desc))
        return results
    except Exception:
        for mol_id, smi in chunk:
            try:
                results.append(predict_single(mol_id, smi))
            except Exception as e:
                results.append({
                    "ID": mol_id, "SMILES": smi, "SIC1": "", "SpMin4_Bhm": "",
                    "SC-5": "", "pIC50": "", "IC50_nM": "", "Activity": "N/A",
                    "status": f"FAILED: {e}",
                })
        return results


def guess_column(fieldnames, candidates):
    lower = {f.lower(): f for f in fieldnames}
    for c in candidates:
        for lf, orig in lower.items():
            if c in lf:
                return orig
    return None


# ----------------------------------------------------------------------
# Visualization helpers
# ----------------------------------------------------------------------

def activity_badge(activity):
    color = ACTIVITY_COLORS.get(activity, "#7f7f7f")
    return f"""
    <span style="
        background-color:{color}22;
        color:{color};
        border:1px solid {color};
        padding:4px 14px;
        border-radius:999px;
        font-weight:600;
        font-size:0.95rem;
    ">{activity}</span>
    """


def render_gauge(pIC50):
    """A gauge showing where the pIC50 sits relative to the activity bands."""
    fig = go.Figure(go.Indicator(
        mode="gauge+number",
        value=pIC50,
        number={"suffix": "  pIC50", "font": {"size": 28}},
        gauge={
            "axis": {"range": [0, 14], "tickwidth": 1},
            "bar": {"color": "#1f2937", "thickness": 0.25},
            "steps": [
                {"range": [0, INACTIVE_CUTOFF], "color": ACTIVITY_COLORS["Inactive"]},
                {"range": [INACTIVE_CUTOFF, MOST_ACTIVE_CUTOFF], "color": ACTIVITY_COLORS["Active"]},
                {"range": [MOST_ACTIVE_CUTOFF, 14], "color": ACTIVITY_COLORS["Most Active"]},
            ],
            "threshold": {
                "line": {"color": "white", "width": 3},
                "thickness": 0.85,
                "value": pIC50,
            },
        },
    ))
    fig.update_layout(height=280, margin=dict(l=20, r=20, t=30, b=10))
    return fig


def render_activity_bar(df):
    counts = (
        df["Activity"].value_counts()
        .reindex(ACTIVITY_ORDER)
        .dropna()
        .reset_index()
    )
    counts.columns = ["Activity", "Count"]
    fig = px.bar(
        counts, x="Activity", y="Count", color="Activity",
        color_discrete_map=ACTIVITY_COLORS, text="Count",
        category_orders={"Activity": ACTIVITY_ORDER},
    )
    fig.update_traces(textposition="outside")
    fig.update_layout(
        showlegend=False, height=380,
        xaxis_title="", yaxis_title="Number of compounds",
        margin=dict(l=10, r=10, t=30, b=10),
    )
    return fig


def render_pic50_distribution(df):
    plot_df = df[pd.to_numeric(df["pIC50"], errors="coerce").notna()].copy()
    plot_df["pIC50"] = plot_df["pIC50"].astype(float)
    fig = px.histogram(
        plot_df, x="pIC50", color="Activity", nbins=30,
        color_discrete_map=ACTIVITY_COLORS,
        category_orders={"Activity": ACTIVITY_ORDER},
        opacity=0.85,
    )
    fig.add_vline(x=INACTIVE_CUTOFF, line_dash="dash", line_color="#d62728",
                   annotation_text="Inactive cutoff (4)", annotation_position="top")
    fig.add_vline(x=MOST_ACTIVE_CUTOFF, line_dash="dash", line_color="#2ca02c",
                   annotation_text="Most active cutoff (11)", annotation_position="top")
    fig.update_layout(
        height=420, bargap=0.05,
        xaxis_title="Predicted pIC50", yaxis_title="Number of compounds",
        legend_title_text="Activity",
        margin=dict(l=10, r=10, t=40, b=10),
    )
    return fig


def render_scatter(df):
    plot_df = df[pd.to_numeric(df["pIC50"], errors="coerce").notna()].copy()
    plot_df["pIC50"] = plot_df["pIC50"].astype(float)
    plot_df["IC50_nM"] = pd.to_numeric(plot_df["IC50_nM"], errors="coerce")
    fig = px.scatter(
        plot_df, x="pIC50", y="IC50_nM", color="Activity",
        color_discrete_map=ACTIVITY_COLORS,
        category_orders={"Activity": ACTIVITY_ORDER},
        hover_data=["ID", "SMILES"], log_y=True,
    )
    fig.update_layout(
        height=420,
        xaxis_title="Predicted pIC50", yaxis_title="IC50 (nM, log scale)",
        legend_title_text="Activity",
        margin=dict(l=10, r=10, t=30, b=10),
    )
    return fig


def render_lipinski_radar(result):
    """Normalized radar chart: each Ro5 property scaled against its own
    threshold, so 1.0 marks the rule-of-five limit line."""
    props = ["MolWt", "LogP", "HBD", "HBA"]
    limits = {"MolWt": 500, "LogP": 5, "HBD": 5, "HBA": 10}
    values = [result[p] / limits[p] for p in props]
    fig = go.Figure()
    fig.add_trace(go.Scatterpolar(
        r=values + [values[0]], theta=props + [props[0]],
        fill="toself", name="Compound",
        line_color="#1f77b4",
    ))
    fig.add_trace(go.Scatterpolar(
        r=[1, 1, 1, 1, 1], theta=props + [props[0]],
        line=dict(color="#d62728", dash="dash"), name="Ro5 limit",
    ))
    fig.update_layout(
        polar=dict(radialaxis=dict(visible=True, range=[0, max(2, max(values) + 0.2)])),
        showlegend=True, height=360,
        margin=dict(l=30, r=30, t=30, b=10),
        title="Rule-of-Five profile (1.0 = threshold)",
    )
    return fig


def render_physchem_scatter(df):
    plot_df = df[pd.to_numeric(df["MolWt"], errors="coerce").notna()].copy()
    plot_df["MolWt"] = plot_df["MolWt"].astype(float)
    plot_df["LogP"] = pd.to_numeric(plot_df["LogP"], errors="coerce")
    fig = px.scatter(
        plot_df, x="MolWt", y="LogP", color="Activity",
        color_discrete_map=ACTIVITY_COLORS,
        category_orders={"Activity": ACTIVITY_ORDER},
        hover_data=["ID", "SMILES", "Lipinski_Violations"],
    )
    fig.add_vline(x=500, line_dash="dash", line_color="gray")
    fig.add_hline(y=5, line_dash="dash", line_color="gray")
    fig.update_layout(
        height=420, xaxis_title="Molecular Weight (Da)", yaxis_title="LogP",
        legend_title_text="Activity", margin=dict(l=10, r=10, t=30, b=10),
    )
    return fig


def style_activity_column(df):
    def _color(val):
        c = ACTIVITY_COLORS.get(val, "#7f7f7f")
        return f"background-color:{c}22; color:{c}; font-weight:600"
    return df.style.applymap(_color, subset=["Activity"])


# ----------------------------------------------------------------------
# Sidebar
# ----------------------------------------------------------------------

with st.sidebar:
    st.markdown("## 🧬 BACE1 QSAR Predictor")
    st.markdown(
        "Predicts pIC50 against **BACE1** (β-secretase 1, an Alzheimer's "
        "disease drug target) from a compound's SMILES, using a GA-MLR "
        "QSAR model built on three PaDEL descriptors: `SIC1`, `SpMin4_Bhm`, "
        "and `SC-5`."
    )
    st.markdown("---")
    st.markdown("### Activity cutoffs")
    st.markdown(
        f"""
        - <span style="color:{ACTIVITY_COLORS['Inactive']}">**Inactive**</span>: pIC50 < {INACTIVE_CUTOFF:g}
        - <span style="color:{ACTIVITY_COLORS['Active']}">**Active**</span>: {INACTIVE_CUTOFF:g} ≤ pIC50 ≤ {MOST_ACTIVE_CUTOFF:g}
        - <span style="color:{ACTIVITY_COLORS['Most Active']}">**Most Active**</span>: pIC50 > {MOST_ACTIVE_CUTOFF:g}
        """,
        unsafe_allow_html=True,
    )
    st.markdown("---")
    mode = st.radio("Mode", ["Single compound", "Batch (CSV)"], index=0)

st.title("BACE1 Inhibitor Activity Predictor")
st.caption(
    "GA-MLR QSAR model · pIC50 → IC50 conversion · activity classification"
)

# ----------------------------------------------------------------------
# Single compound mode
# ----------------------------------------------------------------------

if mode == "Single compound":
    col_in, col_id = st.columns([3, 1])
    with col_in:
        smiles = st.text_input("SMILES", placeholder="e.g. CC(=O)Oc1ccccc1C(=O)O")
    with col_id:
        mol_id = st.text_input("Compound ID (optional)", placeholder="e.g. CMP-01")

    run = st.button("Predict activity", type="primary", use_container_width=False)

    if run:
        if not smiles.strip():
            st.error("Please enter a SMILES string.")
        else:
            with st.spinner("Calculating PaDEL descriptors and scoring..."):
                try:
                    result = predict_single(mol_id.strip() or smiles.strip(), smiles.strip())
                except Exception as e:
                    result = {"status": f"FAILED: {e}"}

            if result["status"].startswith("FAILED"):
                st.error(f"Prediction failed: {result['status']}")
            else:
                if result["status"].startswith("SUSPECT"):
                    st.warning(result["status"])

                st.markdown("### Result")
                c1, c2, c3 = st.columns([1, 1, 1.3])
                with c1:
                    st.metric("Predicted pIC50", f"{result['pIC50']:.3f}")
                with c2:
                    st.metric("IC50 (nM)", f"{result['IC50_nM']:.3g}")
                with c3:
                    st.markdown("**Activity class**")
                    st.markdown(activity_badge(result["Activity"]), unsafe_allow_html=True)

                st.plotly_chart(render_gauge(result["pIC50"]), use_container_width=True)

                if "physchem failed" in result["status"]:
                    st.info(
                        "Physicochemical properties could not be computed for this "
                        "SMILES (RDKit could not parse it)."
                    )
                else:
                    st.markdown("### Physicochemical properties")
                    viol = result["Lipinski_Violations"]
                    p1, p2, p3, p4 = st.columns(4)
                    p1.metric("Molecular Weight", f"{result['MolWt']:.1f} Da")
                    p2.metric("LogP", f"{result['LogP']:.2f}")
                    p3.metric("TPSA", f"{result['TPSA']:.1f} Å²")
                    p4.metric(
                        "Lipinski Ro5 violations", viol,
                        delta="Drug-like" if viol == 0 else f"{viol} broken",
                        delta_color="normal" if viol == 0 else "inverse",
                    )
                    st.plotly_chart(render_lipinski_radar(result), use_container_width=True)
                    with st.expander("All physicochemical values"):
                        st.dataframe(
                            pd.DataFrame([{c: result[c] for c in PHYSCHEM_COLUMNS}]),
                            use_container_width=True, hide_index=True,
                        )

                with st.expander("QSAR descriptor values"):
                    st.dataframe(
                        pd.DataFrame([{
                            "ID": result["ID"], "SMILES": result["SMILES"],
                            "SIC1": result["SIC1"], "SpMin4_Bhm": result["SpMin4_Bhm"],
                            "SC-5": result["SC-5"],
                        }]),
                        use_container_width=True, hide_index=True,
                    )

# ----------------------------------------------------------------------
# Batch mode
# ----------------------------------------------------------------------

else:
    st.markdown(
        "Upload a CSV with a **SMILES** column (and optionally an **ID** column). "
        "The app will guess the columns automatically, or you can pick them below."
    )
    uploaded = st.file_uploader("CSV file", type=["csv"])

    if uploaded is not None:
        try:
            text = uploaded.getvalue().decode("utf-8-sig")
            reader = csv.DictReader(io.StringIO(text))
            fieldnames = reader.fieldnames or []
            if not fieldnames:
                st.error("Could not read a header row from this CSV.")
                st.stop()
            all_rows = list(reader)
        except UnicodeDecodeError:
            st.error("Could not decode the file as UTF-8 — please upload a plain CSV.")
            st.stop()
        except csv.Error as e:
            st.error(f"Malformed CSV: {e}")
            st.stop()

        guessed_id = guess_column(fieldnames, ["id", "sid", "name"])
        guessed_smiles = guess_column(fieldnames, ["smile"])

        c1, c2 = st.columns(2)
        with c1:
            id_col = st.selectbox(
                "ID column", fieldnames,
                index=fieldnames.index(guessed_id) if guessed_id in fieldnames else 0,
            )
        with c2:
            smiles_options = fieldnames
            default_idx = smiles_options.index(guessed_smiles) if guessed_smiles in smiles_options else 0
            smiles_col = st.selectbox("SMILES column", smiles_options, index=default_idx)

        rows = []
        for row in all_rows:
            smi = (row.get(smiles_col) or "").strip()
            mid = (row.get(id_col) or "").strip() if id_col else smi
            if smi:
                rows.append((mid or smi, smi))

        st.caption(f"{len(rows)} compound(s) with a non-empty SMILES value found.")

        if st.button("Run batch prediction", type="primary", disabled=(len(rows) == 0)):
            progress = st.progress(0, text="Starting...")
            all_results = []
            total = len(rows)
            try:
                for start in range(0, total, CHUNK_SIZE):
                    chunk = rows[start:start + CHUNK_SIZE]
                    all_results.extend(process_chunk(chunk))
                    done = len(all_results)
                    progress.progress(
                        min(done / total, 1.0),
                        text=f"Processing {done} / {total} molecules...",
                    )
                progress.progress(1.0, text="Done.")
            except Exception as e:
                traceback.print_exc()
                st.error(f"Batch run failed: {e}")
                st.stop()

            df = pd.DataFrame(all_results)
            st.session_state["batch_df"] = df

    if "batch_df" in st.session_state:
        df = st.session_state["batch_df"]

        n_total = len(df)
        n_ok = (df["status"] == "ok").sum()
        n_suspect = df["status"].astype(str).str.startswith("SUSPECT").sum()
        n_fail = df["status"].astype(str).str.startswith("FAILED").sum()

        st.markdown("### Summary")
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Total compounds", n_total)
        m2.metric("Scored OK", int(n_ok))
        m3.metric("Suspect", int(n_suspect))
        m4.metric("Failed", int(n_fail))

        act_counts = df["Activity"].value_counts()
        a1, a2, a3 = st.columns(3)
        a1.metric("🔴 Inactive", int(act_counts.get("Inactive", 0)))
        a2.metric("🟡 Active", int(act_counts.get("Active", 0)))
        a3.metric("🟢 Most Active", int(act_counts.get("Most Active", 0)))

        n_druglike = (pd.to_numeric(df["Lipinski_Violations"], errors="coerce") == 0).sum()
        st.caption(f"💊 {int(n_druglike)} / {n_total} compounds are fully Lipinski Ro5-compliant (0 violations).")

        st.markdown("### Activity — QSAR results")
        v1, v2 = st.columns(2)
        with v1:
            st.plotly_chart(render_activity_bar(df), use_container_width=True)
        with v2:
            st.plotly_chart(render_pic50_distribution(df), use_container_width=True)

        st.plotly_chart(render_scatter(df), use_container_width=True)

        st.markdown("### Physicochemical properties (Lipinski Ro5)")
        st.plotly_chart(render_physchem_scatter(df), use_container_width=True)

        st.markdown("### Results table")
        filter_choice = st.multiselect(
            "Filter by activity", ACTIVITY_ORDER, default=ACTIVITY_ORDER
        )
        shown = df[df["Activity"].isin(filter_choice)] if filter_choice else df
        st.dataframe(style_activity_column(shown), use_container_width=True, hide_index=True)

        csv_bytes = df.to_csv(index=False).encode("utf-8")
        st.download_button(
            "⬇️ Download full results (CSV)",
            data=csv_bytes,
            file_name="BACE1_predictions.csv",
            mime="text/csv",
            type="primary",
        )
