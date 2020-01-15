# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'pyEPR'
copyright = '2020, Zlatko Minev, Zaki Leghtas, and the pyEPR Team'
author = 'Zlatko Minev, Zaki Leghtas, and the pyEPR Team'

import sys
import os
sys.path.insert(0, os.path.abspath("../../pyEPR"))
print(sys.path)

# The full version, including alpha/beta/rc tags
import pyEPR
version = pyEPR.__version__
release = version

import sphinx_rtd_theme

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
	"sphinx.ext.intersphinx",
	"sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    'sphinx.ext.coverage', 
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    "sphinx.ext.mathjax",
    "sphinx_rtd_theme",
    "IPython.sphinxext.ipython_directive",
    "IPython.sphinxext.ipython_console_highlighting",
    "matplotlib.sphinxext.plot_directive",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["**.ipynb_checkpoints"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme =  'sphinx_rtd_theme' #'default' # 'sphinx_rtd_theme' #'alabaster' "sphinxdoc" 'classic'
if 0:
	import os
	on_rtd = os.environ.get('READTHEDOCS') == 'True'
	if on_rtd:
	    html_theme = 'default'
	else:
	    html_theme = 'nature'


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


# If false, no module index is generated.
html_use_modindex = True

html_show_sourcelink = True

# Sort members by type
autodoc_member_order = 'groupwise'

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"