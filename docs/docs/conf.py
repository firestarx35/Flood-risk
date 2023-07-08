import sys
import os


# -- Project information -----------------------------------------------------

project = "Team Ouse's Flood Tool"
copyright = 'Team Ouse'
author = 'Team Ouse'

# The full version, including alpha/beta/rc tags
release = '1.0'

sys.path.insert(0, os.path.abspath(os.sep.join((os.curdir,'..'))))

project = 'Flood Tool'
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinx.ext.viewcode']
source_suffix = '.rst'
master_doc = 'index'
exclude_patterns = ['_build']
autoclass_content = "both"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'