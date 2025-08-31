# src/settings.py
from __future__ import annotations
import copy
import streamlit as st

# -------- Defaults για όλα τα βήματα --------
DEFAULTS = {
    "preprocessing": {
        "min_counts_cells": 1000,
        "min_cells_genes": 3,
        "target_sum": 1e4,
        "n_top_genes": 2000,
        "scale_max": 10,
        "n_pcs": 40,
        "n_neighbors": 15,
    },
    "integration": {
        "method": "scanorama",
        "n_pcs": 40,
        "n_neighbors": 15,
    },
    "annotation": {
        "method": "leiden",
        "resolution": 1.0,
    },
    "deg": {
        "method": "wilcoxon",
        "padj_threshold": 0.05,
        "log2fc_threshold": 1.0,
    },
    "plots": {
        "umap_point_size": 8,
        "volcano_annotate": 0,  # 0: χωρίς labels, >0: top N
    },
}

# -------- State helpers --------
def ensure_settings():
    if "settings" not in st.session_state:
        st.session_state["settings"] = copy.deepcopy(DEFAULTS)

def reset_settings():
    st.session_state["settings"] = copy.deepcopy(DEFAULTS)

# -------- Sidebar UI --------
def settings_sidebar(sections: tuple[str, ...]):
    """Sidebar ρυθμίσεις, μόνο για τα ζητούμενα sections."""
    ensure_settings()
    s = st.session_state["settings"]

    if "preprocessing" in sections:
        with st.sidebar.expander("🧹 Preprocessing", expanded=False):
            pp = s["preprocessing"]
            pp["min_counts_cells"] = st.number_input("Min counts / cell", 0, 1_000_000, pp["min_counts_cells"])
            pp["min_cells_genes"]  = st.number_input("Min cells / gene", 0, 50_000,   pp["min_cells_genes"])
            pp["target_sum"]       = st.number_input("Normalize target_sum", 0.0, 1e8, pp["target_sum"], step=1e3, format="%.0f")
            pp["n_top_genes"]      = st.number_input("HVGs (n_top_genes)", 100, 20_000, pp["n_top_genes"], step=100)
            pp["scale_max"]        = st.number_input("Scale max value", 1, 100, pp["scale_max"])
            c1, c2 = st.columns(2)
            pp["n_pcs"]       = c1.number_input("PCA components", 2, 200, pp["n_pcs"])
            pp["n_neighbors"] = c2.number_input("Neighbors (graph)", 2, 200, pp["n_neighbors"])

    if "integration" in sections:
        with st.sidebar.expander("🔗 Integration", expanded=False):
            it = s["integration"]
            it["method"]      = st.selectbox("Method", ["scanorama"], index=0)
            c1, c2 = st.columns(2)
            it["n_pcs"]       = c1.number_input("PCA components", 2, 200, it["n_pcs"])
            it["n_neighbors"] = c2.number_input("Neighbors (graph)", 2, 200, it["n_neighbors"])

    if "annotation" in sections:
        with st.sidebar.expander("🏷️ Annotation", expanded=False):
            an = s["annotation"]
            an["method"]     = st.selectbox("Method", ["leiden"], index=0)
            an["resolution"] = st.slider("Leiden resolution", 0.1, 3.0, an["resolution"], 0.1)

    if "deg" in sections:
        with st.sidebar.expander("📊 DEG", expanded=False):
            dg = s["deg"]
            dg["method"]           = st.selectbox("Method", ["wilcoxon", "t-test", "t-test_overestim_var", "logreg"],
                                                   index=["wilcoxon","t-test","t-test_overestim_var","logreg"].index(dg["method"]))
            c1, c2 = st.columns(2)
            dg["padj_threshold"]   = c1.slider("padj threshold", 0.0, 0.5, dg["padj_threshold"], 0.01)
            dg["log2fc_threshold"] = c2.slider("|log2FC| threshold", 0.0, 5.0, dg["log2fc_threshold"], 0.1)

    if "plots" in sections:
        with st.sidebar.expander("📈 Plots", expanded=False):
            pl = s["plots"]
            pl["umap_point_size"]  = st.slider("UMAP point size", 2, 20, pl["umap_point_size"])
            pl["volcano_annotate"] = st.number_input("Volcano annotate top N", 0, 100, pl["volcano_annotate"])

    with st.sidebar:
        if st.button("↺ Επαναφορά προεπιλογών"):
            reset_settings()
            st.rerun()

    return s

# -------- Summary chips (προαιρετικό visual) --------
def _chip(label: str, value) -> str:
    return f"<span style='display:inline-block;border:1px solid rgba(0,0,0,.12);padding:4px 8px;border-radius:999px;margin:2px 6px 0 0;font-size:12px'>{label}: <b>{value}</b></span>"

def render_summary(section_keys: tuple[str, ...]):
    ensure_settings()
    s = st.session_state["settings"]
    html = "<div style='margin: 6px 0 10px 0'>"
    if "preprocessing" in section_keys:
        pp = s["preprocessing"]
        html += _chip("HVGs", pp["n_top_genes"]) + _chip("PCs", pp["n_pcs"]) + _chip("Neighbors", pp["n_neighbors"])
    if "integration" in section_keys:
        it = s["integration"]
        html += _chip("Integration", it["method"]) + _chip("PCs", it["n_pcs"])
    if "annotation" in section_keys:
        an = s["annotation"]
        html += _chip("Annot", an["method"]) + _chip("Res", an["resolution"])
    if "deg" in section_keys:
        dg = s["deg"]
        html += _chip("DEG", dg["method"]) + _chip("padj", dg["padj_threshold"]) + _chip("|log2FC|", dg["log2fc_threshold"])
    if "plots" in section_keys:
        pl = s["plots"]
        html += _chip("UMAP size", pl["umap_point_size"]) + _chip("Volcano N", pl["volcano_annotate"])
    html += "</div>"
    st.markdown(html, unsafe_allow_html=True)
