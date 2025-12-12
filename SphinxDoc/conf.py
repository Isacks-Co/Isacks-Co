# Configuration file for the Sphinx documentation builder.

import os
import sys
from datetime import datetime

# ------------------------------------------------------------------
# Path setup: make your project source importable
# ------------------------------------------------------------------
sys.path.insert(0, os.path.abspath("../SourceCode"))

# Optional extra paths
# sys.path.insert(0, os.path.abspath("../Scripts"))
# sys.path.insert(0, os.path.abspath("../DefectPlots"))

# Mock heavy / external imports so autodoc doesn't crash
autodoc_mock_imports = [
    "torch",
    "mace",
    "httk",
    "ase",
    "numpy",
    "pandas",
    "asap3",
    "classes",
    "matplotlib",
    "scipy",
    "abad_classes",
    "pymongo",
    "optimade",
    "DataBase_scripts"
]

# ------------------------------------------------------------------
# Project information
# ------------------------------------------------------------------
project = "Computational physics project - TFYA99"
author = "Isacks & Co"
current_year = datetime.now().year
copyright = f"{current_year}, {author}"
release = "0.1.0"

# ------------------------------------------------------------------
# General configuration
# ------------------------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

autosummary_generate = True

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True

autodoc_typehints = "description"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

pygments_style = "solarized-light"

# ------------------------------------------------------------------
# HTML output
# ------------------------------------------------------------------
html_theme = "sphinx_rtd_theme"

# ------------------------------------------------------------------
# LaTeX / PDF output (plain, no boxes)
# ------------------------------------------------------------------
latex_engine = "pdflatex"

latex_elements = {
    "papersize": "a4paper",
    "pointsize": "11pt",

    # Clean, standard font
    "fontpkg": r"""
\usepackage{lmodern}
""",

    # Minimal preamble
    "preamble": r"""
\usepackage{microtype}
\usepackage{xcolor}

% Simple readable links
\hypersetup{
    colorlinks=true,
    linkcolor=blue,
    urlcolor=blue,
    citecolor=blue,
}

% Better verbatim wrapping
\usepackage{fvextra}
\DefineVerbatimEnvironment{Highlighting}{Verbatim}{
    breaklines,
    breakanywhere,
    commandchars=\\\{\},
    fontsize=\small
}
""",
}

latex_documents = [
    (
        "index",
        "SphinxDocumentation.tex",
        "Computational physics project - TFYA99",
        "Isacks & Co",
        "manual",
    ),
]

latex_domain_indices = True
