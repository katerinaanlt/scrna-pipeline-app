"""Entry point: ``streamlit run main.py``."""
import streamlit as st

from scrna_app.pages import PAGES, sidebar_summary

st.set_page_config(page_title="scRNA-seq Pipeline", page_icon="🧬", layout="wide")
st.sidebar.title("🧬 scRNA-seq Pipeline")
choice = st.sidebar.radio("Navigation", list(PAGES))
sidebar_summary()
PAGES[choice]()
