import streamlit as st
from src import utils

st.set_page_config(page_title="Expression", page_icon="🎚️", layout="wide")
st.title("🎚️ Expression plots")

adata = (st.session_state.get("adata_annot")
         or st.session_state.get("adata_integrated")
         or st.session_state.get("adata_merged"))
if adata is None:
    st.info("Χρειάζεσαι δεδομένα από προηγούμενες σελίδες.")
    st.stop()

gene = st.text_input("Γονίδιο", value="")
groupby = st.text_input("Ομαδοποίηση (π.χ. 'cell_type' ή άδειο)", value="cell_type")

if st.button("Σχεδίαση"):
    if not gene.strip():
        st.warning("Γράψε ένα γονίδιο.")
    else:
        try:
            fig = utils.violin_expression_figure(adata, gene.strip(), groupby.strip() or None)
            st.pyplot(fig, width="stretch")
        except Exception as e:
            st.error(f"Σφάλμα: {e}")
