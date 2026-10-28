"""
BACE1 QSAR Visualizations
-------------------------
Charts and display helpers used by bace1_qsar_app.py.
"""

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

# Activity classification thresholds (on predicted pIC50)
INACTIVE_CUTOFF = 4.0
MOST_ACTIVE_CUTOFF = 11.0

ACTIVITY_COLORS = {
    "Inactive": "#d62728",
    "Active": "#f2c744",
    "Most Active": "#2ca02c",
    "N/A": "#7f7f7f",
}
ACTIVITY_ORDER = ["Inactive", "Active", "Most Active", "N/A"]


def activity_badge(activity):
    """Return an HTML badge coloured by activity class."""
    color = ACTIVITY_COLORS.get(activity, "#7f7f7f")
    style = (
        f"background-color:{color}22; color:{color}; "
        f"border:1px solid {color}; padding:4px 14px; "
        "border-radius:999px; font-weight:600; font-size:0.95rem;"
    )
    return f"<span style='{style}'>{activity}</span>"


def render_gauge(pIC50):
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
    """Each Ro5 property scaled to its own limit, so 1.0 = the threshold."""
    props = ["MolWt", "LogP", "HBD", "HBA"]
    limits = {"MolWt": 500, "LogP": 5, "HBD": 5, "HBA": 10}
    values = [result[p] / limits[p] for p in props]
    fig = go.Figure()
    fig.add_trace(go.Scatterpolar(
        r=values + [values[0]], theta=props + [props[0]],
        fill="toself", name="Compound", line_color="#1f77b4",
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

    styler = df.style
    # Styler.map needs pandas >= 2.1; older versions use applymap
    apply = styler.map if hasattr(styler, "map") else styler.applymap
    return apply(_color, subset=["Activity"])