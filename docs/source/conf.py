# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath("../../"))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'PoseButcher'
copyright = '2024, Max Winokan'
author = 'Max Winokan'
github_username = 'mwinokan'
github_repository = 'PoseButcher'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
              'sphinx.ext.autodoc',
              'autoapi.extension',
              'sphinx.ext.doctest',
              'sphinx.ext.extlinks',
              'sphinx.ext.mathjax',
              'sphinx.ext.viewcode',
              'sphinx.ext.napoleon',
              'sphinx.ext.intersphinx']

autoapi_dirs = ['../../posebutcher/']

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

source_suffix = '.rst'
master_doc = 'index'
# project = 'ASE'
# copyright = f'{datetime.date.today().year}, ASE-developers'
# templates_path = ['templates']
exclude_patterns = ['build']
# default_role = 'math'
# pygments_style = 'sphinx'
autoclass_content = 'both'
modindex_common_prefix = ['mp.']
