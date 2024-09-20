# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
import os
import sys
import datetime
from importlib import import_module


# -- Project information -----------------------------------------------------

project = "Ramses"
author = "Romain Teyssier"
copyright = '{0}, {1}'.format(datetime.datetime.now().year, "CEA and Romain Teyssier")
package = "Ramses"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "myst_parser",
    "sphinx.ext.mathjax",
    "nbsphinx",
    "sphinx_copybutton",
    'sphinx_simplepdf',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints", 
                    "wiki/Content.md", "README.md", "wiki/Home.md"]

source_suffix = ['.rst', '.md']
# source_suffix = ".rst"
html_sourcelink_suffix = ""  # Avoid .ipynb.txt extensions in sources

# The master toctree document.
master_doc = "index"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
primary_color = '#000000'
primary_color_rgb_string = 'rgba(0, 0, 0, 1)' # For pdf, same color as `primary_color`
secondary_color = '#FFD587'

html_theme = "sphinx_book_theme"
html_theme_options = {
    "repository_url": "https://github.com/ramses-organisation/ramses",
    "repository_branch": "main",
    "use_repository_button": True,
    "use_issues_button": True,
    "use_edit_page_button": True,
    "show_toc_level": 3,
}
html_logo = "./img/logo.svg"
html_favicon = "./img/logo.svg"

html_title = '{0}'.format(project)

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

nbsphinx_timeout = 300


suppress_warnings = []


simplepdf_vars = {
    'primary': primary_color,
    'secondary': secondary_color,
    'cover': '#ffffff',
    'white': '#ffffff',
    'links': 'FA2323',
    'cover-bg': primary_color_rgb_string,
    'top-left-content': 'counter(page)',
}
