# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath("../.."))


project = 'FInDS'
copyright = '2026, The University of Houston'
author = 'Mufaro J. Machaya'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",  # optional but recommended
    "sphinx.ext.imgmath",
    "sphinxcontrib.bibtex"
]

# latex config
latex_documents = [
    (
        'index',                 # 1. Start doc (usually 'index')
        'finds.tex',             # 2. Target name (This determines the PDF filename!)
        'FInDS: A Fish Interaction and Dynamics Simulator', # 3. Document Title
        'Mufaro J. Machaya',     # 4. Author
        'manual'                 # 5. Document class (e.g., manual, howto)
    ),
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

napoleon_use_ivar = True

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
