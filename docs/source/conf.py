# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath("../.."))

import src.main
dir(src.main)

project = 'fish'
copyright = '2026, Mufaro J. Machaya'
author = 'Mufaro J. Machaya'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",  # optional but recommended
    "sphinx.ext.imgmath",
    "sphinxcontrib.bibtex"
]

bibtex_bibfiles = ["references.bib"]
bibtex_default_style = 'plain'

templates_path = ['_templates']
exclude_patterns = []

#autodoc configuration
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "member-order": "bysource",
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
