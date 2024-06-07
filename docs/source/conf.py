# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import os
import sys
sys.path.insert(0, os.path.abspath('/Users/benjaminyuen/Documents/Documents - Benjamin’s MacBook Pro/Ben/physics/Birmingham/QuantumCode/diracpy_project/src/diracpy'))

# Path to the Graphviz dot executable
graphviz_dot = 'dot'  # Assuming 'dot' is in your PATH

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Diracpy'
copyright = '2024, Ben Yuen'
author = 'Ben Yuen'
release = '1.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',  # If you use Google or NumPy style docstrings
    'sphinx.ext.viewcode',  # To add links to highlighted source code
    'sphinx.ext.inheritance_diagram',
    # 'sphinx.ext.graphviz',
    'sphinxcontrib.bibtex',
    'sphinx_rtd_theme',
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'
# html_theme = 'furo'
# html_theme = 'piccolo_theme'
html_permalinks_icon = '§'
# html_theme = 'insipid'
html_static_path = ['_static']

# -- Extension configuration -------------------------------------------------
# Configure the bibtex_bibfiles setting
bibtex_bibfiles = ['references.bib']
