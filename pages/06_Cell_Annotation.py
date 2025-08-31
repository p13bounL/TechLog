import streamlit as st
from src import utils
from src.settings import ensure_settings, settings_sidebar, render_summary

st.set_page_config(page_title="Cell Annotation", page_icon="🏷️", layout="wide")
st.title("🏷️ Cell Annotation")

ensure_settings()
settings = settings_sidebar(sections=("annotation","plots"))
render_summary(("annotation","plots"))

base = (st.session_state.get("adata_integrated") or st.session_state.get("adata_merged"))
if base is None:
    st.info("Χρειάζεσαι δεδομένα από 04 ή/και 05.")
    st.stop()

if st.button("Τρέξε annotation / επανεκτέλεση"):
    with st.spinner("Τρέχει annotation..."):
        ann = utils.annotate_cells(base, resolution=settings["annotation"]["resolution"])
        st.session_state["adata_annot"] = ann
    st.success("Ολοκληρώθηκε το annotation.")

ann = st.session_state.get("adata_annot")
if ann is not None:
    if "cell_type" in ann.obs:
        st.write("Κατανομή cell_type:")
        st.write(ann.obs["cell_type"].astype(str).value_counts().head(20))

    cols = list(ann.obs.columns)
    if not cols:
        st.warning("Δεν υπάρχουν στήλες στο `.obs` για χρωματισμό.")
    else:
        default_val = "cell_type" if "cell_type" in cols else cols[0]
        default_index = cols.index(default_val)
        color_by = st.selectbox(
            "Χρωματισμός UMAP ανά:",
            options=cols,
            index=default_index if "annot_umap_color" not in st.session_state else
                  (cols.index(st.session_state["annot_umap_color"]) if st.session_state["annot_umap_color"] in cols else default_index),
            key="annot_umap_color",
        )
        ann_plot = utils.ensure_umap(ann)  # μπορείς να περάσεις n_pcs/n_neighbors αν θέλεις να είναι ίδιες με integration
        fig = utils.umap_figure(ann_plot, color=color_by, point_size=settings["plots"]["umap_point_size"])
        st.pyplot(fig, width="stretch")
else:
    st.info("Πάτα το κουμπί πιο πάνω για να τρέξει το annotation.")
