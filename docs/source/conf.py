# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Sphinx Access to project code -------------------------------------------

import os
import sys
# Добавляем путь к проекту в sys.path
sys.path.insert(0, os.path.abspath('../commented_code'))
# Проверяем правильность пути
print("Path to project:", os.path.abspath('../commented_code'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Bioformats Dream'
copyright = '2025, Tsarkov Alexander, Shabalin Artyom, Lipinskaya Alyona, Shkuratova Valentina'
author = 'Tsarkov Alexander, Shabalin Artyom, Lipinskaya Alyona, Shkuratova Valentina'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinxcontrib.mermaid'
]

templates_path = ['_templates']
exclude_patterns = []

language = 'ru'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
