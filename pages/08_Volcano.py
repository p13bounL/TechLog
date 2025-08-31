import streamlit as st
from src import utils
from src.settings import ensure_settings, settings_sidebar, render_summary

st.set_page_config(page_title="Volcano", page_icon="🌋", layout="wide")
st.title("🌋 Volcano Plot")

ensure_settings()
settings = settings_sidebar(sections=("deg","plots"))
render_summary(("deg","plots"))

df = st.session_state.get("deg_df")
if df is None:
    st.info("Δεν υπάρχει DEG πίνακας. Πήγαινε στη σελίδα 08.")
    st.stop()

fig = utils.volcano_figure(
    df,
    log2fc_col="log2FC",
    p_col="padj",
    fc_thr=settings["deg"]["log2fc_threshold"],
    p_thr=settings["deg"]["padj_threshold"],
    annotate_top=settings["plots"]["volcano_annotate"],
)
st.pyplot(fig, width="stretch")
