import streamlit as st
from textwrap import dedent

st.set_page_config(page_title="Members", page_icon="👥", layout="wide")

# ----------------- DATA: Συμπληρώνουμε εδώ τα μέλη της ομάδας -----------------
MEMBERS = [
    {
        "name": "Λάμπρος Μπουντούλης",
        "role": "Επί Πτυχίω - Π2013075",
        "email": "p13boun@ionio.gr",
        "github": "https://github.com/username",
        "photo": "",  # π.χ. "https://avatars.githubusercontent.com/u/9919"
    },
    # + πρoσθήκη κι άλλων μελών...
]

# ----------------- STYLES -----------------
st.markdown(dedent("""
<style>
.cards { display: grid; grid-template-columns: repeat( auto-fit, minmax(260px, 1fr) ); gap: 16px;  margin-top: 12px; }
.card {
  background: var(--background-color);
  border: 1px solid rgba(0,0,0,.08);
  border-radius: 16px;
  padding: 16px;
  box-shadow: 0 4px 16px rgba(0,0,0,.04);
  transition: transform .15s ease, box-shadow .15s ease;
}
.card:hover { transform: translateY(-2px); box-shadow: 0 8px 24px rgba(0,0,0,.06); }
.header { display:flex; align-items:center; gap:12px; }
.avatar {
  width:56px; height:56px; border-radius:50%;
  display:flex; align-items:center; justify-content:center;
  background: linear-gradient(135deg, #e0f2fe, #f1f5f9);
  color:#0f172a; font-weight:700; font-size:18px;
  border: 1px solid rgba(0,0,0,.08);
}
.name   { font-size:18px; font-weight:700; line-height:1.1; margin:0; }
.role   { font-size:13px; opacity:.8; margin:2px 0 0 0; }
.about  { font-size:13px; margin:10px 0 0 0; }
.links  { display:flex; gap:10px; margin-top:10px; flex-wrap:wrap; }
.link a { text-decoration:none; font-size:13px; border:1px solid rgba(0,0,0,.12); padding:6px 8px; border-radius:999px; }
.link a:hover { background: rgba(0,0,0,.04); }
.badge { display:inline-flex; align-items:center; gap:6px; }
.badge .icon { font-size:14px; }
.section-title { margin-top: 6px; }
</style>
"""), unsafe_allow_html=True)

# ----------------- HEADER -----------------
st.title("👥 Members")
st.caption("Η ομάδα που υλοποίησε την εφαρμογή.")

# ----------------- FILTERS (προαιρετικό)
cols = st.columns(3)
with cols[0]:
    q = st.text_input("🔎 Αναζήτηση", value="")
with cols[1]:
    show_email = st.toggle("Εμφάνιση email", value=True)
with cols[2]:
    sort_by = st.selectbox("Ταξινόμηση κατά", ["Σειρά λίστας", "Όνομα", "Ρόλο"], index=0)

# ----------------- FILTER / SORT LOGIC -----------------
items = MEMBERS[:]
if q.strip():
    ql = q.lower()
    def hit(m):
        return any(ql in str(m.get(k,"")).lower() for k in ["name","role","about"])
    items = [m for m in items if hit(m)]

if sort_by == "Όνομα":
    items = sorted(items, key=lambda m: m.get("name","").lower())
elif sort_by == "Ρόλο":
    items = sorted(items, key=lambda m: m.get("role","").lower())

# ----------------- RENDER -----------------
def initials(name: str) -> str:
    parts = [p for p in name.split() if p]
    if not parts: return "👤"
    if len(parts) == 1: return parts[0][:2].upper()
    return (parts[0][0] + parts[-1][0]).upper()

st.markdown('<div class="cards">', unsafe_allow_html=True)

for m in items:
    name = m.get("name","")
    role = m.get("role","")
    email = m.get("email","")
    gh = m.get("github","")
    li = m.get("linkedin","")
    photo = m.get("photo","").strip()
    about = m.get("about","")

    # avatar: είτε εικόνα, είτε αρχικά
    if photo:
        avatar_html = f'<img src="{photo}" class="avatar" alt="avatar" />'
    else:
        avatar_html = f'<div class="avatar">{initials(name)}</div>'

    links_html = ""
    links = []
    if gh:
        links.append(f'<span class="link"><a href="{gh}" target="_blank" rel="noopener">🐙 GitHub</a></span>')
    if li:
        links.append(f'<span class="link"><a href="{li}" target="_blank" rel="noopener">🔗 LinkedIn</a></span>')
    if show_email and email:
        links.append(f'<span class="link"><a href="mailto:{email}">✉️ {email}</a></span>')
    if links:
        links_html = f'<div class="links">{"".join(links)}</div>'

    card = f"""
    <div class="card">
      <div class="header">
        {avatar_html}
        <div>
          <p class="name">{name}</p>
          <p class="role">{role}</p>
        </div>
      </div>
      <p class="about">{about}</p>
      {links_html}
    </div>
    """
    st.markdown(card, unsafe_allow_html=True)

st.markdown("</div>", unsafe_allow_html=True)

# ----------------- FOOTNOTE -----------------
st.divider()
st.caption("Πρόσθεσε/άλλαξε μέλη στη λίστα `MEMBERS` στην κορυφή του αρχείου.")
