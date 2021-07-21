# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import datetime
import sphinx_rtd_theme

# -- Project information -----------------------------------------------------
year = datetime.date.today().year
project = "MANDYOC"
copyright = "2020-{}, The {} Developers".format(year, project)

# The full version, including alpha/beta/rc tags
release = "0.1.0"


# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.graphviz",
    "sphinxcontrib.bibtex",
    # 'sphinx.ext.autosectionlabel',
    # 	'sphinx.ext.autodoc',
    # 	'sphinx.ext.intersphinx',
    # 	'sphinx-prompt',
    # 	'recommonmark',
    # 	'notfound.extension',
    # 	'hoverxref.extension',
    # 	'sphinx_search.extension',
    # 	'sphinxemoji.sphinxemoji',
]

graphviz_output_format = "svg"

bibtex_bibfiles = ["files/refs.bib"]
bibtex_default_style = "unsrt"
# bibtex_reference_style = 'author_year'

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]
# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix = ".rst"
master_doc = "index"

# -- Options for HTML output -------------------------------------------------
# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_theme_options = {
    "logo_only": False,
    "display_version": True,
    "prev_next_buttons_location": "bottom",
    "vcs_pageview_mode": "",
    # 'style_nav_header_background': 'white',
    # Toc options
    "collapse_navigation": True,
    "sticky_navigation": True,
    "navigation_depth": 4,
    "includehidden": True,
    "titles_only": False,
}

html_last_updated_fmt = "%b %d, %Y"
html_title = project
html_short_title = project
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
pygments_style = "sphinx"


def setup(app):
    app.add_css_file("css/custom.css")


numfig = True
math_numfig = True
numfig_secnum_depth = 1
math_eqref_format = "Eq. {number}"
