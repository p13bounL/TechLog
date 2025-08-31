import streamlit as st
from src import utils

st.set_page_config(page_title="Preprocessing", page_icon="🧹", layout="wide")
st.title("🧹 Προεπεξεργασία .h5ad (batch filtering)")

st.markdown("Δώσε φακέλους όπως στο notebook σου (π.χ. `./data/h5ad` → `./data/h5ad_filt`).")

colA, colB = st.columns(2)
with colA:
    input_dir = st.text_input("Φάκελος εισόδου (.h5ad)", value="./data/h5ad")
with colB:
    output_dir = st.text_input("Φάκελος εξόδου", value="./data/h5ad_filt")

col1, col2, col3 = st.columns(3)
with col1:
    n_min = st.number_input("n_genes_min", min_value=0, value=1000, step=100)
with col2:
    n_max = st.number_input("n_genes_max", min_value=1, value=10000, step=100)
with col3:
    compression = st.selectbox("compression", ["zstd", "gzip", "lzf"], index=0)

if st.button("Τρέξε batch preprocessing"):
    try:
        written = utils.batch_preprocess_h5ad_dir(
            input_dir=input_dir,
            output_dir=output_dir,
            n_genes_min=int(n_min),
            n_genes_max=int(n_max),
            compression=compression,
        )
        if written:
            st.success(f"Αποθηκεύτηκαν {len(written)} αρχεία.")
            st.write(written)
        else:
            st.info("Δεν βρέθηκαν .h5ad στον φάκελο εισόδου.")
    except Exception as e:
        st.error(f"Σφάλμα: {e}")
