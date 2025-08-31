import streamlit as st
from src import utils
from src.settings import ensure_settings, settings_sidebar, render_summary

st.set_page_config(page_title="Integration", page_icon="🔗", layout="wide")
st.title("🔗 Integration")

ensure_settings()
settings = settings_sidebar(sections=("integration","plots"))
render_summary(("integration","plots"))

merged = st.session_state.get("adata_merged")
if merged is None:
    st.info("Δεν υπάρχει merged AnnData. Πήγαινε στη σελίδα 04 για ένωση.")
    st.stop()

if st.button("Τρέξε integration / επανεκτέλεση"):
    with st.spinner("Τρέχει integration..."):
        # εδώ μπορείς να περάσεις τις παραμέτρους στο δικό σου integrate_scanorama αν το υποστηρίξει
        integ = utils.integrate_scanorama(merged)
        st.session_state["adata_integrated"] = integ
    st.success("Ολοκληρώθηκε το integration.")

integ = st.session_state.get("adata_integrated")
if integ is not None:
    cols = list(integ.obs.columns)
    if not cols:
        st.warning("Δεν υπάρχουν στήλες στο `.obs` για χρωματισμό.")
    else:
        default_val = "batch" if "batch" in cols else cols[0]
        default_index = cols.index(default_val)
        # ΧΩΡΙΣ programmatic set στο session_state πριν το widget (φεύγει το warning)
        color_by = st.selectbox(
            "Χρωματισμός UMAP ανά:",
            options=cols,
            index=default_index if "integ_umap_color" not in st.session_state else
                  (cols.index(st.session_state["integ_umap_color"]) if st.session_state["integ_umap_color"] in cols else default_index),
            key="integ_umap_color",
        )
        # UMAP με παραμέτρους από settings
        integ_plot = utils.ensure_umap(
            integ,
            n_neighbors=settings["integration"]["n_neighbors"],
            n_pcs=settings["integration"]["n_pcs"],
        )
        fig = utils.umap_figure(integ_plot, color=color_by, point_size=settings["plots"]["umap_point_size"])
        st.pyplot(fig, width="stretch")
else:
    st.info("Πάτα το κουμπί πιο πάνω για να τρέξει το integration.")
