import streamlit as st
from src import utils

st.set_page_config(page_title="Intro", page_icon="📘", layout="wide")
st.title("📘 Intro - Έλεγχος Περιβάλλοντος")

# Έλεγχος εκδόσεων (προστασία αν λείπει η συνάρτηση)
try:
    info = utils.check_versions()
    st.success("Έλεγχος εκδόσεων OK")
    st.json(info)
except Exception as e:
    st.info(f"Παράκαμψη check_versions (utils): {e}")

# Δίσκος
try:
    free_gb = utils.disk_free()
    st.markdown(f"**Ελεύθερος χώρος δίσκου:** ~{free_gb:.1f} GB")
except Exception:
    pass

# Εξερεύνηση φακέλων
st.subheader("📂 Εξερεύνηση φακέλων/αρχείων")
path = st.text_input("Δώσε φάκελο ή αρχείο για να δεις τι υπάρχει", value=".")
if st.button("Λίστα περιεχομένων"):
    try:
        items = utils.list_dir(path)
        st.write(items)
    except Exception as e:
        st.error(f"Σφάλμα: {e}")
