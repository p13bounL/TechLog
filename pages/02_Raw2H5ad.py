from pathlib import Path
import streamlit as st
from src import utils

st.set_page_config(page_title="Raw → H5AD", page_icon="📦", layout="wide")
st.title("📦 Μετατροπή από 10x mtx φάκελο σε .h5ad")

# ---- Βοηθητικό dropdown για να μην γράφεις path με το χέρι ----
base = Path("./data/RAW")
cands = [str(p) for p in base.iterdir() if p.is_dir()] if base.exists() else []
pick = st.selectbox(
    "Διαθέσιμοι φάκελοι (./data/RAW)",
    ["(χειροκίνητα)"] + cands,
    index=0,
    key="p02_pick",
)

folder = st.text_input(
    "Ή γράψε πλήρη διαδρομή (μέσα στο container)",
    value="" if pick == "(χειροκίνητα)" else pick,
    key="p02_folder",
)

out_path = st.text_input(
    "Όνομα αρχείου εξόδου (.h5ad)",
    value="./data/h5ad/sample01.h5ad",
    key="p02_outpath",
)

condition = st.text_input(
    "Προαιρετικό condition (π.χ. control/disease)",
    value="",
    key="p02_condition",
)

col1, col2 = st.columns(2)
with col1:
    if st.button("Διάβασε φάκελο", key="p02_read_btn"):
        if folder.strip():
            try:
                adata = utils.read_10x_mtx_folder(
                    folder.strip(),
                    var_names="gene_symbols",
                    condition=(condition or None),
                )
                st.session_state["adata_raw"] = adata
                st.success(f"Φορτώθηκε! σχήμα: {adata.shape}")
                st.write("Obs columns:", list(adata.obs.columns)[:10])
                st.write("Var columns:", list(adata.var.columns)[:10])
            except Exception as e:
                st.error(str(e))
        else:
            st.error("Γράψε διαδρομή φακέλου.")

with col2:
    if st.button("Αποθήκευση σε .h5ad", key="p02_save_btn"):
        adata = st.session_state.get("adata_raw")
        if adata is None:
            st.warning("Πρώτα φόρτωσε τα δεδομένα (αριστερό κουμπί).")
        else:
            path = utils.save_h5ad(adata, out_path.strip())
            st.success(f"Αποθηκεύτηκε: {path}")
