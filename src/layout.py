import streamlit as st

def sidebar_members_bottom_button(
    label: str = "Members",
    icon: str = "👥",
    key: str = "members_btn_app",
) -> bool:
    """Κουμπί στο ΚΑΤΩ μέρος του sidebar (φαίνεται μόνο όταν είναι open)."""
    st.markdown("""
    <style>
    [data-testid="stSidebar"] > div:first-child { display:flex; flex-direction:column; height:100%; }
    .sb-bottom-wrap { position: sticky; bottom: 8px; padding: 8px; }
    :root:not(:has([data-testid="stSidebar"][aria-expanded="true"])) .sb-bottom-wrap { display:none !important; }
    .sb-bottom-btn { width:100%; padding:10px 12px; border-radius:12px;
                     border:1px solid rgba(0,0,0,.12); background:transparent; }
    .sb-bottom-btn:hover { background: rgba(0,0,0,.06); }
    </style>
    """, unsafe_allow_html=True)

    with st.sidebar:
        st.markdown('<div class="sb-bottom-wrap">', unsafe_allow_html=True)
        clicked = st.button(f"{icon} {label}", key=key)
        st.markdown('</div>', unsafe_allow_html=True)
    return clicked
