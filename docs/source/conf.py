# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import sys, os
import subprocess

sys.path.insert(0, os.path.abspath('../../src'))

project = 'MOXελ'
copyright = '2023-2024, Antonios P. Sarikas'
author = 'Antonios P. Sarikas'
release = f'{subprocess.run(["git", "describe"])}'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
        'sphinx.ext.autodoc',
        'sphinx.ext.viewcode',
        'sphinx.ext.napoleon',
        'sphinxemoji.sphinxemoji',
        'sphinx_copybutton',
        'sphinxarg.ext',
        'sphinx_code_tabs',
        ]

templates_path = ['_templates']
exclude_patterns = ['modules.rst']

# Exclude input prompts from copybutton
copybutton_exclude = '.linenos, .gp, .go'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
#html_static_path = ['_static']
