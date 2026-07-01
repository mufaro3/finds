# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath("../.."))

project = 'FINDS'
copyright = '2026, The University of Houston'
author = 'Mufaro J. Machaya'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",  # optional but recommended
    "sphinxcontrib.bibtex"
]

# latex config
latex_documents = [
    (
        'index',                 # 1. Start doc (usually 'index')
        'finds.tex',             # 2. Target name (This determines the PDF filename!)
        'FINDS: The Fish INteraction and Dynamics Simulator', # 3. Document Title
        'Mufaro J. Machaya',     # 4. Author
        'manual'                 # 5. Document class (e.g., manual, howto)
    ),
]

latex_elements = {
    "extraclassoptions": "openany",
    "pointsize": "12pt",
    'passoptionstopackages': r'''
\PassOptionsToPackage{svgnames}{xcolor}
''',
    'fontpkg': r'''
\usepackage{tgpagella}
\usepackage{tgheros}
\usepackage{inconsolata}
''',
    'preamble': r'''
\usepackage[titles]{tocloft}
\cftsetpnumwidth {1.25cm}\cftsetrmarg{1.5cm}
\setlength{\cftchapnumwidth}{0.75cm}
\setlength{\cftsecindent}{\cftchapnumwidth}
\setlength{\cftsecnumwidth}{1.25cm}
''',
    'fncychap': r'\usepackage[Bjornstrup]{fncychap}',
    'printindex': r'\footnotesize\raggedright\printindex',
}

bibtex_bibfiles = ["references.bib"]
bibtex_default_style = 'plain'

templates_path = ['_templates']
exclude_patterns = []

#autodoc configuration
autodoc_typehints = "none"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "member-order": "bysource",
}
autoclass_content = "class"
napoleon_use_ivar = True

# signature folding
HIDE_SIGNATURES = {
    "finds.io.IO",
    "finds.postprocessing.process_data",
    "finds.processingmodules.DensityAnimationGenerator.DensityAnimationGenerator",
    "finds.processingmodules.AnimationGenerator.AnimationGenerator",
    "finds.processingmodules.CrossSectionGenerator.CrossSectionGenerator",
    "finds.processingmodules.MeanRadialDistancePlot.MeanRadialDistancePlot"
}

def process_signature(
    app,
    what,
    name,
    obj,
    options,
    signature,
    return_annotation,
):
    if name in HIDE_SIGNATURES:
        return "(...)", return_annotation

    return signature, return_annotation


def setup(app):
    app.connect(
        "autodoc-process-signature",
        process_signature
    )

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
html_favicon = '_static/favicon.ico'
