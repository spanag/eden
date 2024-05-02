# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os, sys
from pathlib import Path

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

repo_local_path = '..'
repo_local_path = os.path.abspath(repo_local_path)

project = 'EDEN'
copyright = ' 2019– , EDEN authors and contributors'
author = 'EDEN authors and contributors'

import eden_simulator
release = eden_simulator.__version__
# release = Path(repo_local_path+"/VERSION").read_text().strip()
print(release)

# get a version tag or branch, for cross references to the git repo
ver_tag = 'v'+release
if ver_tag in os.listdir(repo_local_path+'/.git/refs/tags'): repo_tag = ver_tag
else: repo_tag = 'development'

# https://stackoverflow.com/questions/27381997/get-variables-in-sphinx-templates
html_context = {}

repo_gitlab_user = 'c7859/neurocomputing-lab/Inferior_OliveEMC'
repo_github_user = 'spanag'
repo_name = 'eden'
repo_url_prefix = f'https://gitlab.com/{repo_gitlab_user}/{repo_name}'
repo_url_prefix_github = f'https://github.com/{repo_github_user}/{repo_name}'
repo_gitlab_user_binder = repo_gitlab_user.replace('/','%2F')
line_length_limit = 100_000_000_000

for x in ['repo_gitlab_user', 'repo_github_user', 'repo_gitlab_user_binder', 'repo_name', 'repo_url_prefix', 'repo_url_prefix_github', 'repo_tag','line_length_limit']: html_context[x] = globals()[x]

print(repo_local_path) 
print(html_context)
html_theme_optionss = {'lalala':'444'}

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# add the python package(s) to python path

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
# print(f"path is {os.environ['EDEN_PYTHON_PACKAGE_PATH']}")
# sys.path.insert(0, os.path.abspath(os.environ['EDEN_PYTHON_PACKAGE_PATH']))


extensions = [
	'sphinx.ext.autodoc',
	'sphinx.ext.autosummary',
	'sphinx.ext.napoleon',
	
    # "sphinx.ext.doctest",
	# 'sphinx_gallery.load_style',
	'nbsphinx',
	"sphinx_design",
	
]

templates_path = ['_templates']
exclude_patterns = ['Thumbs.db', '.DS_Store']

# autodoc
# Mock unavailable library modules
autodoc_mock_imports = ["lxml", "numpy",] # could import these for now, but not necessarily in the general case...
# autodoc_default_flags = ['members', 'undoc-members' ]

# Document __init__, __repr__, and __str__ methods? why not 
# https://www.woolseyworkshop.com/2020/07/17/documenting-python-programs-with-sphinx/
def autodoc_skip_member(app, what, name, obj, would_skip, options):
	if name in ("__init__", "__repr__", "__str__"):
		return False
	return would_skip

autodoc_default_options = {
    "exclude-members": "main, parse_dict_arg, parse_list_arg, build_namespace, convert_case, process_args",
    "imported-members":True,
}

autosummary_generate = True  # Turn on sphinx.ext.autosummary, because why have it just work https://stackoverflow.com/questions/62613202/automatically-document-all-modules-recursively-with-sphinx-autodoc
# it needs autyomodule to not be used on the same page, to show the module (if automodule precedes it then autosummary applies to children of the module)
autosummary_ignore___all__ = False # sane behaviour was added in late 2021, as a niche alternative. But it doens't have an effect anyway...
	
# nbsphinx

# hardcoded thumbnails for notebooks which normally don't show them, and for .gif thumbnails
nbsphinx_thumbnails = {
    'tut_net': '_static/thumb_tut_net.gif',
    'exa_lfp': '_static/thumb_exa_lfp.png',
}

# This is processed by Jinja2 and inserted before each notebook https://github.com/spatialaudio/nbsphinx/blob/0.9.3/doc/conf.py#L43
# https://github.com/spatialaudio/nbsphinx/issues/419#issuecomment-603522609
# Add links to Binder, TODO deep, colab
# TODO labpath instead of filepath?
# TODO .binder !
# TODO https://www.sphinx-doc.org/en/master/usage/theming.html#builtin-themes
# navigation_with_keys  body_max_width https://pydata-sphinx-theme.readthedocs.io/en/stable/user_guide/layout.html#references

nbsphinx_prolog = r"""
{% set docname = 'docs/' + env.doc2path(env.docname, base=None) %}
{% set v = env.config.html_context %}

.. raw:: html

    <div id="run-notebook-online" class="admonition tip">
      This page was generated from
      <a class="reference external" href="{{ v.repo_url_prefix|e }}/blob/{{ v.repo_tag|e }}/{{ docname|e }}">{{ docname|e }}</a>.
      <span style="white-space: nowrap;">Interactive online version:</span>
      
	  <span style="white-space: nowrap;"><a href="https://mybinder.org/v2/gl/{{ v.repo_gitlab_user_binder|e }}%2F{{ v.repo_name }}/{{ v.repo_tag|e }}?filepath={{ docname|e }}"><img alt="Binder badge" src="https://mybinder.org/badge_logo.svg" style="vertical-align:text-bottom"></a></span>
	  
	  <span style="white-space: nowrap;"><a href="https://colab.research.google.com/github/{{ v.repo_github_user }}/{{ v.repo_name }}/blob/{{ v.repo_tag|e }}/{{ docname|e }}"><img alt="Colab badge" src="https://colab.research.google.com/assets/colab-badge.svg" style="vertical-align:text-bottom"></a></span>
	  
	  <span style="white-space: nowrap;"><a href="https://deepnote.com/launch?url={{ v.repo_url_prefix_github|e }}/blob/{{ v.repo_tag|e }}/{{ docname|e }}"><img alt="Deepnote badge" src="https://img.shields.io/badge/launch-deepnote-3793EF?logo=Deepnote" style="vertical-align:text-bottom"></a>.</span>
	  
      <span style="white-space: nowrap;"><a href="{{ env.docname.split('/')|last|e + '.ipynb' }}" class="reference download internal" download>Download notebook</a>.</span>
    </div>

.. raw:: latex

    \nbsphinxstartnotebook{\scriptsize\noindent\strut
    \textcolor{gray}{The following section was generated from
    \sphinxcode{\sphinxupquote{\strut {{ docname | escape_latex }}}} \dotfill}}
"""

# This is processed by Jinja2 and inserted after each notebook
nbsphinx_epilog = r"""
{% set docname = 'docs/' + env.doc2path(env.docname, base=None) %}
.. raw:: latex

    \nbsphinxstopnotebook{\scriptsize\noindent\strut
    \textcolor{gray}{\dotfill\ \sphinxcode{\sphinxupquote{\strut
    {{ docname | escape_latex }}}} ends here.}}
"""

# common
rst_prolog = r"""
.. |contact us| replace:: :ref:`contact us <contact_us>`
.. |Contact us| replace:: :ref:`Contact us <contact_us>`

.. https://stackoverflow.com/questions/19686897/how-can-i-link-to-a-page-section-in-a-sphinx-toctree
.. role:: hidden
    :class: hidden

"""

def setup(app):
	import os
	os.makedirs("docs/_images", exist_ok=True)
	
	# autodoc
	app.connect("autodoc-skip-member", autodoc_skip_member)
	def builder_inited(app):
		env = app.env
		env.settings['line_length_limit'] = 10_000_000_000 # NB: there's not much use making files bigger than the default 100M though
	
	app.connect('builder-inited', builder_inited)

# -- Options for EPUB output --

# These are just defined to avoid Sphinx warnings related to EPUB: (from nbsphinx repo)
version = release
suppress_warnings = ['epub.unknown_project_files']

# -- Options for HTML output --
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static'] # TODO add gathered image files?
# html_css_files = ['_static/custom.css'] # HOTE: this does not seem to work, the css file is not even copied. But at least custom.css works for alabaster and pydata theme

# TODO html_logo, favicon
# TODO get from static store
# TODO html_sidebars and nml link ! https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-html_sidebars

