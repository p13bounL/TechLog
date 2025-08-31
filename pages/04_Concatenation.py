import streamlit as st
import scanpy as sc
import tempfile
from src import utils

st.set_page_config(page_title="Concatenation", page_icon="🧩", layout="wide")
st.title("🧩 Συνένωση .h5ad (Concatenation)")

files = st.file_uploader("Φόρτωσε ΠΟΛΛΑ .h5ad", type=["h5ad"], accept_multiple_files=True)
save_path = st.text_input("(Προαιρετικό) Αρχείο εξόδου για merged (.h5ad)", value="")

if st.button("Ένωση αρχείων"):
    if not files:
        st.warning("Διάλεξε αρχεία πρώτα.")
    else:
        try:
            adatas = []
            for uf in files:
                with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
                    tmp.write(uf.getbuffer())
                    tmp_path = tmp.name
                adatas.append(sc.read_h5ad(tmp_path))

            merged = utils.concatenate_adatas(adatas)
            st.session_state["adata_merged"] = merged
            st.success(f"OK! merged shape: {merged.shape}")
            if save_path.strip():
                utils.save_h5ad(merged, save_path.strip())
                st.info(f"Αποθηκεύτηκε: {save_path.strip()}")
        except Exception as e:
            st.error(f"Σφάλμα: {e}")
