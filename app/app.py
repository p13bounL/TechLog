import streamlit as st
import time
from src.layout import sidebar_members_bottom_button
from src.members_ui import render_members

st.set_page_config(page_title="App", page_icon="🧬", layout="wide")

# ---------- Members dialog ----------
@st.dialog("👥 Members")
def _members_dialog():
    render_members()

# ---------- Sidebar Members button (πάντα στην αρχή, πριν οτιδήποτε άλλο) ----------
if sidebar_members_bottom_button():
    _members_dialog()

# ---------- CSS: γωνιακά emojis + shift όταν ανοίγει το sidebar ----------
SYMBOL_SIZE_PX = 72
SYMBOL_OPACITY = 0.16
TOP_OFFSET = 16
SIDE_OFFSET = 16
SIDEBAR_OPEN_PX = 300

st.markdown(f"""
<style>
:root {{
  --corner-size:{SYMBOL_SIZE_PX}px; --corner-opacity:{SYMBOL_OPACITY};
  --top-off:{TOP_OFFSET}px; --side-off:{SIDE_OFFSET}px; --sbw:0px; --sbw-open:{SIDEBAR_OPEN_PX}px;
}}
:root:has([data-testid="stSidebar"][aria-expanded="true"]) {{ --sbw:var(--sbw-open); }}

.corner {{ position:fixed; font-size:var(--corner-size); opacity:var(--corner-opacity);
  pointer-events:none; z-index:0; filter:drop-shadow(0 2px 4px rgba(0,0,0,.25)); }}
.corner.tl {{ top:var(--top-off); left:calc(var(--side-off) + var(--sbw)); animation: float 6s ease-in-out infinite; }}
.corner.bl {{ bottom:var(--top-off); left:calc(var(--side-off) + var(--sbw)); animation: float 8s ease-in-out infinite; }}
.corner.tr {{ top:var(--top-off); right:var(--side-off); animation: float 7s ease-in-out infinite; }}
.corner.br {{ bottom:var(--top-off); right:var(--side-off); animation: float 9s ease-in-out infinite; }}
@keyframes float {{ 0%{{transform:translateY(0) rotate(0)}} 50%{{transform:translateY(-6px) rotate(2deg)}} 100%{{transform:translateY(0) rotate(0)}} }}
.hero {{ font-size:26px; line-height:1.6; text-align:center; margin-top:14vh; z-index:1; position:relative; }}
.cursor {{ display:inline-block; width:1ch; animation:blink 1s steps(2,start) infinite; }}
@keyframes blink {{ to {{ visibility:hidden; }} }}
</style>
<div class="corner tl">🧬</div><div class="corner tr">🔬</div>
<div class="corner bl">🧫</div><div class="corner br">🧪</div>
""", unsafe_allow_html=True)

# ---------- Typing effect: μόνο την πρώτη φορά ----------
TEXT = 'Εργασία Εξαμήνου για το Μάθημα "Τεχνολογία Λογισμικού". Ανάπτυξη εφαρμογής για εκτέλεση διεργασιών ανάλυσης δεδομένων και Μηχανικής Μάθησης σε δεδομένα Μοριακής Βιολογίας.'
if "hero_done" not in st.session_state:
    st.session_state["hero_done"] = False

if not st.session_state["hero_done"]:
    speed = st.slider("Ταχύτητα πληκτρολόγησης (χαρακτ./δευτ.)", 10, 60, 30, 1)
    delay = 1.0 / speed
    ph = st.empty()
    typed = ""
    for ch in TEXT:
        typed += ch
        ph.markdown(f"<div class='hero'>{typed}<span class='cursor'>|</span></div>", unsafe_allow_html=True)
        time.sleep(delay)
    st.session_state["hero_done"] = True
else:
    st.markdown(f"<div class='hero'>{TEXT}</div>", unsafe_allow_html=True)

if st.button("↺ Ξαναπαίξε"):
    st.session_state["hero_done"] = False
    st.rerun()
