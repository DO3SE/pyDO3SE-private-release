"""Entrypoint for pyDO3SE GUI."""
import streamlit as st

from pyDO3SE_gui.plugin_examples.photosynthesis.ewert import (
    gui_ewert_co2_concentration_in_stomata_iteration,
    gui_ewert, gui_ewert_co2_concentration_in_stomata_loop,
    gui_ewert_full_year,
)


def home():
    """Home page."""
    st.title("Home")
    st.write('Pick a page from the left hand side')


if __name__ == "__main__":
    pages = {
        "home": home,
        "ewert_co2_concentration_in_stomata_iteration":
        gui_ewert_co2_concentration_in_stomata_iteration,
        "ewert_co2_concentration_in_stomata_loop":
        gui_ewert_co2_concentration_in_stomata_loop,
        "ewert_single_hour": gui_ewert,
        "ewert_full_year": gui_ewert_full_year,
    }
    page = st.sidebar.selectbox("Choose Page", list(pages.keys()))
    pages[page]()
