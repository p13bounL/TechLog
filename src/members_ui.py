import streamlit as st
from textwrap import dedent
import base64, mimetypes
from pathlib import Path

# >>> Συμπλήρωσε τα δικά σας στοιχεία εδώ <<<
MEMBERS = [
    {"name":"Λάμπρος Μπουντούλης","role":"Επί Πτυχίω - Π2013075","email":"p13boun@ionio.gr","github":"","photo":"",},
    # ...
]

def _photo_src(photo_str: str | None) -> str | None:
    if not photo_str: return None
    ps = photo_str.strip()
    if ps.startswith(("http://","https://","data:")):
        return ps
    p = Path(ps)
    if not p.exists():
        p = Path.cwd() / ps
        if not p.exists():
            return None
    mime = mimetypes.guess_type(p.name)[0] or "image/png"
    b64 = base64.b64encode(p.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{b64}"

def render_members():
    st.markdown(dedent("""
    <style>
    .cards { display:grid; grid-template-columns:repeat(auto-fit,minmax(260px,1fr)); gap:16px; margin-top:8px; }
    .card { border:1px solid rgba(0,0,0,.08); border-radius:16px; padding:16px; box-shadow:0 4px 16px rgba(0,0,0,.04); }
    .header { display:flex; align-items:center; gap:12px; }
    .avatar { width:56px; height:56px; border-radius:50%; display:flex; align-items:center; justify-content:center;
              background:linear-gradient(135deg,#e0f2fe,#f1f5f9); font-weight:700; }
    .name{margin:0;font-weight:700} .role{margin:2px 0 0 0;opacity:.8;font-size:13px}
    .links{display:flex;gap:10px;margin-top:10px;flex-wrap:wrap}
    .links a{border:1px solid rgba(0,0,0,.12);padding:6px 8px;border-radius:999px;text-decoration:none;font-size:13px}
    </style>
    """), unsafe_allow_html=True)

    st.subheader("👥 Members")
    st.caption("Η ομάδα που υλοποίησε την εφαρμογή.")

    def initials(name: str) -> str:
        parts = [p for p in name.split() if p]
        if not parts: return "👤"
        return (parts[0][0] + (parts[-1][0] if len(parts)>1 else "")).upper()

    st.markdown('<div class="cards">', unsafe_allow_html=True)
    for m in MEMBERS:
        src = _photo_src((m.get("photo") or "").strip())
        avatar = f'<img src="{src}" class="avatar" alt="avatar"/>' if src \
                 else f'<div class="avatar">{initials(m.get("name",""))}</div>'
        links=[]
        if m.get("github"):   links.append(f'<a href="{m["github"]}" target="_blank">🐙 GitHub</a>')
        if m.get("linkedin"): links.append(f'<a href="{m["linkedin"]}" target="_blank">🔗 LinkedIn</a>')
        if m.get("email"):    links.append(f'<a href="mailto:{m["email"]}">✉️ {m["email"]}</a>')
        links_html = " ".join(links)
        st.markdown(
            f"""<div class="card">
                   <div class="header">{avatar}
                     <div><p class="name">{m.get("name","")}</p>
                          <p class="role">{m.get("role","")}</p></div>
                   </div>
                   <p style="margin-top:10px">{m.get("about","")}</p>
                   <div class="links">{links_html}</div>
                </div>""",
            unsafe_allow_html=True,
        )
    st.markdown('</div>', unsafe_allow_html=True)
