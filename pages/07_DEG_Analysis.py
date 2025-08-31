import streamlit as st
import numpy as np
import pandas as pd
from src import utils
from src.settings import ensure_settings, settings_sidebar, render_summary

st.set_page_config(page_title="DEG", page_icon="📊", layout="wide")
st.title("📊 Differential Expression (2 ομάδες)")

ensure_settings()
settings = settings_sidebar(sections=("deg",))
render_summary(("deg",))

# Πηγή δεδομένων
adata_src = (st.session_state.get("adata_annot")
             or st.session_state.get("adata_integrated")
             or st.session_state.get("adata_merged"))
if adata_src is None:
    st.info("Χρειάζεσαι δεδομένα από προηγούμενες σελίδες.")
    st.stop()

# 1) Επιλογή στήλης groupby
obs_cols = list(adata_src.obs.columns)
groupby = st.selectbox(
    "Στήλη groupby",
    options=obs_cols,
    index=obs_cols.index("cell_type") if "cell_type" in obs_cols else 0
)

col = adata_src.obs[groupby]

# Θα χρησιμοποιήσουμε ένα αντίγραφο αν χρειαστεί να δημιουργήσουμε προσωρινή στήλη
adata_deg = adata_src
groupby_used = groupby
group1 = None
group2 = None

# 2) Αν η στήλη είναι συνεχής (numeric) → φτιάξε ΔΥΟ ομάδες με split
if pd.api.types.is_numeric_dtype(col):
    st.info("Η επιλεγμένη στήλη είναι συνεχής. Διάλεξε τον τρόπο διαχωρισμού σε 2 ομάδες.")
    mode = st.radio("Τύπος διαχωρισμού", ["Median split", "Quantiles (Q1 vs Q4)", "Manual threshold"], horizontal=True)

    tmp_col = "__tmp_deg_group__"
    adata_deg = adata_src.copy()
    values = adata_deg.obs[groupby].astype(float)

    if mode == "Median split":
        med = float(np.nanmedian(values))
        left_label  = f"≤ {med:.3g}"
        right_label = f"> {med:.3g}"
        adata_deg.obs[tmp_col] = np.where(values <= med, left_label, right_label)
        group1, group2 = left_label, right_label

    elif mode == "Quantiles (Q1 vs Q4)":
        q = st.slider("Ποσοστό για άκρα", 0.05, 0.45, 0.25, 0.05)
        q1 = float(np.nanquantile(values, q))
        q4 = float(np.nanquantile(values, 1 - q))
        low_label  = f"<= Q{int(q*100)} ({q1:.3g})"
        high_label = f">= Q{int((1-q)*100)} ({q4:.3g})"
        group = np.full(values.shape[0], None, dtype=object)
        group[values <= q1] = low_label
        group[values >= q4] = high_label
        adata_deg.obs[tmp_col] = group
        # πέτα το “μεσαίο” 50% ώστε να μείνουν καθαρές οι δύο άκρες
        adata_deg = adata_deg[pd.notna(adata_deg.obs[tmp_col])].copy()
        group1, group2 = low_label, high_label

    else:  # Manual threshold
        mn, mx = float(np.nanmin(values)), float(np.nanmax(values))
        default_thr = float(np.nanmedian(values))
        thr = st.number_input("Κατώφλι (threshold)", mn, mx, default_thr)
        left_label  = f"≤ {thr:.3g}"
        right_label = f"> {thr:.3g}"
        adata_deg.obs[tmp_col] = np.where(values <= thr, left_label, right_label)
        group1, group2 = left_label, right_label

    groupby_used = tmp_col

# 3) Αν η στήλη είναι κατηγορική/διακριτή → διάλεξε δύο κατηγορίες
else:
    if pd.api.types.is_categorical_dtype(col):
        options = list(col.cat.categories)
    else:
        options = sorted(pd.unique(col.dropna()))
    group1 = st.selectbox("Ομάδα A", options=options, format_func=str, key="deg_group1")
    group2 = st.selectbox("Ομάδα B", options=[o for o in options if o != group1] or options, format_func=str, key="deg_group2")

# 4) Μέγεθος ομάδων (πριν τρέξουμε)
n1 = int((adata_deg.obs[groupby_used] == group1).sum())
n2 = int((adata_deg.obs[groupby_used] == group2).sum())
st.caption(f"Μέγεθος ομάδων: A={n1} | B={n2}")

# 5) Τρέξε DEG
if st.button("Τρέξε DEG"):
    if group1 == group2:
        st.error("Οι δύο ομάδες πρέπει να είναι διαφορετικές.")
        st.stop()
    if n1 < 2 or n2 < 2:
        st.error(f"Πολύ μικρές ομάδες για DEG: A={n1}, B={n2}. Διάλεξε άλλη στήλη ή διαφορετικό split.")
        st.stop()

    try:
        df = utils.deg_two_groups(
            adata_deg,
            groupby_used,
            group1,
            group2,
            method=settings["deg"]["method"],
        )
        st.session_state["deg_df"] = df
        st.success(f"ΟΚ, {len(df)} γραμμές.")
        st.dataframe(df.head(100), width="stretch")
        if not df.empty:
            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button("Λήψη ως CSV", data=csv, file_name="deg_results.csv", mime="text/csv")
    except Exception as e:
        st.error(f"DEG error: {e}")
