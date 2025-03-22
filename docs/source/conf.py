# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import sys, os
from importlib.metadata import version as get_version

sys.path.insert(0, os.path.abspath('../../src'))

project = 'MOXελ'
copyright = '2023-2024, Antonios P. Sarikas'
author = 'Antonios P. Sarikas'
release = get_version('pymoxel')

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
        'sphinx.ext.autodoc',
        'sphinx.ext.viewcode',
        'sphinx.ext.napoleon',
        'sphinx.ext.intersphinx',
        'sphinx_copybutton',
        'sphinx_design',
        'sphinx_issues',
        ]

#templates_path = ['_templates']
#exclude_patterns = ['modules.rst']
autodoc_typehints = 'none'

# Exclude input prompts from copybutton
copybutton_exclude = '.linenos, .gp, .go'

# The package is too small.
#intersphinx_mapping = {
#        'python': ('https://docs.python.org/3', None),
#        'numpy': ('https://numpy.org/doc/stable/', None),
#        }

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_logo = 'images/moxel_logo.svg'
html_theme_options = {'logo_only': True, 'version_selector': False}
html_static_path = ['_static']
html_css_files = ['custom.css']

# Path to GitHub repo {group}/{project}  (note that `group` is the GitHub user or organization)
issues_github_path = "adosar/moxel"

# which is the equivalent to:
issues_uri = "https://github.com/{group}/{project}/issues/{issue}"
issues_prefix = "#"
issues_pr_uri = "https://github.com/{group}/{project}/pull/{pr}"
issues_pr_prefix = "#"
issues_commit_uri = "https://github.com/{group}/{project}/commit/{commit}"
issues_commit_prefix = "@"
issues_user_uri = "https://github.com/{user}"
issues_user_prefix = "@"
