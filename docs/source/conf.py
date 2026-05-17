# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path information -----------------------------------------------------

import sys
import os

sys.path.insert(0, os.path.abspath("../../pyEPR"))
print(sys.path)

# -- Project information -----------------------------------------------------

project = "pyEPR"
copyright = "2017-2025, Zlatko Minev, Zaki Leghtas, and the pyEPR Team"
author = "Zlatko Minev, Zaki Leghtas, and the pyEPR Team"

# html_title controls the browser tab and search-engine title.
# We omit the version number — it dates quickly and clutters search results.
html_title = "Welcome to pyEPR! — Energy-Participation-Ratio Framework"

# The full version, including alpha/beta/rc tags
import pyEPR

version = pyEPR.__version__
release = version

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.intersphinx",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",  #  parse both NumPy and Google style docstrings
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    "sphinx.ext.mathjax",
    "sphinx_design",
    "myst_nb",
    "matplotlib.sphinxext.plot_directive",
]

nb_execution_mode = "off"   # use saved outputs, never re-execute
myst_enable_extensions = ["colon_fence", "dollarmath"]

# https://github.com/readthedocs/readthedocs.org/issues/2569
master_doc = "index"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["**.ipynb_checkpoints"]

suppress_warnings = [
    "myst.header",               # non-consecutive header levels in old notebooks
    "misc.highlighting_failure", # IPython ? magic in notebook code cells
    "ref.python",                # duplicate QuantumAnalysis cross-reference in reports.rst
]


# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = False

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
show_authors = True

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"


numpydoc_show_class_members = True
napoleon_numpy_docstring = True
napoleon_use_admonition_for_notes = True

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "pydata_sphinx_theme"

html_theme_options = {
    "github_url": "https://github.com/zlatko-minev/pyEPR",
    "use_edit_page_button": False,
    "show_toc_level": 2,
    "navbar_align": "left",
    "navbar_end": ["navbar-icon-links"],
    "secondary_sidebar_items": ["page-toc"],
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version"],
    "icon_links": [
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/pyEPR-quantum/",
            "icon": "fa-solid fa-box",
        },
    ],
}
# Add any paths that contain custom themes here, relative to this directory.


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]


# If false, no module index is generated.
html_use_modindex = True

html_show_sourcelink = True

# Sort members by type
# autodoc_member_order = 'groupwise'

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"


# -----------------------------------------------------------------------------
# Autosummary
# -----------------------------------------------------------------------------

autosummary_generate = True

# -----------------------------------------------------------------------------
# Autodoc
# -----------------------------------------------------------------------------
# The supported options are
# 'members', 'member-order', 'undoc-members', 'private-members',
# 'special-members', 'inherited-members', 'show-inheritance', 'ignore-module-all',
# 'imported-members' and 'exclude-members'.
autodoc_default_options = {
    "inherited-members": None,
    #'member-order': 'bysource',
    "member-order": "alphabetical",  # This value selects if automatically documented members are sorted alphabetical (value 'alphabetical'), by member type (value 'groupwise') or by source order (value 'bysource'). The default is alphabetical.
    "undoc-members": True,  # Members without docstrings will be left out, unless you give the undoc-members flag option:
    "exclude-members": "__weakref__",
    "show-inheritance": True,  # , a list of base classes will be inserted just below the class signature (when used with automodule, this will be inserted for every class that is documented in the module).
}


# If true, figures, tables and code-blocks are automatically numbered if they
# have a caption.
numfig = True

# A dictionary mapping 'figure', 'table', 'code-block' and 'section' to
# strings that are used for format of figure numbers. As a special character,
# %s will be replaced to figure number.
numfig_format = {"table": "Table %s"}
# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = "en"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "colorful"

# A boolean that decides whether module names are prepended to all object names
# (for object types where a “module” of some kind is defined), e.g. for
# py:function directives.
add_module_names = True

# A list of prefixes that are ignored for sorting the Python module index
# (e.g., if this is set to ['foo.'], then foo.bar is shown under B, not F).
# This can be handy if you document a project that consists of a single
# package. Works only for the HTML builder currently.
# modindex_common_prefix = ['pyEPR.']
