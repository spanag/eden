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

import os
confdir = os.getcwd()

# get binary assets from static store like html_logo, favicon, LATER put somewhere appropriate
def get_more_assets(app):
	import urllib3
	eden_assets_prefix = 'https://eden-simulator.org/assets/docs/'
	eden_assets_files = [
		'eden_logo_white_bg.png',
		'thumb_intro_neuroml.png',
		'thumb_intro_lems.png',
		'thumb_extension_pointers.png',
		'thumb_extension_io.png',
		'thumb_extension_writable.png',
		'thumb_extension_multiflux.png',
		'favicon.png',
	]
	if app.builder.format == 'latex': # LATER check if latex...
		eden_assets_files += [
			'tutorial_network_balls.png',
			'tutorial_network_neuron.png',
			'tutorial_network_detailed.png',
			'example_lfp_3d_vm.png',
			'example_lfp_3d_current.png',
			'example_lfp_3d_full.png',
			'extension_customsetup_balls.png',
		]

	target_prefix = confdir+'/_static/'
	urls_files = [ (eden_assets_prefix+x, target_prefix+x) for x in eden_assets_files ]

	http = urllib3.PoolManager()
	# https://www.owenrumney.co.uk/retry-urllib3-requests/
	retry = urllib3.util.Retry(3, redirect=20, raise_on_status=True, status_forcelist=range(400, 600))
	for url, filename in urls_files:
		print(url)
		data = http.request('GET', url, retries=retry, timeout=60).data
		# print('%r page is %d bytes' % (url, len(data)))
		with open(filename, "wb") as f: f.write(data) 

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
    'sphinx_codeautolink',  # automatic links from code to documentation NEXT
	
	'sphinx.ext.intersphinx',
	# "sphinx.ext.autosectionlabel", # more trouble than it's worth, for multiple chapters
	
	# "sphinx.ext.doctest",
	# 'sphinx_gallery.load_style',
	'nbsphinx',
	'sphinx_design',
	'sphinx_copybutton',
	
    # 'sphinxcontrib.bibtex',  # for bibliographic references?
    'sphinxcontrib.rsvgconverter',  # for SVG->PDF conversion in LaTeX output
]

templates_path = ['_templates']
exclude_patterns = ['Thumbs.db', '.DS_Store']

suppress_warnings = [
	'codeautolink.match_name', # old python? probably not? https://github.com/felix-hilden/sphinx-codeautolink/issues/96
	'docutils', # https://github.com/spatialaudio/nbsphinx/issues/753
]

# intersphinx
intersphinx_mapping = {
	'neuroml': ('https://docs.neuroml.org/', None),
	'neuron': ('https://nrn.readthedocs.io/en/latest/', None),
	'k3d': ('https://k3d-jupyter.org/', None),
	'pyvista': ('https://docs.pyvista.org', None),
}

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

add_module_names = False # at least for print...

# Source code links
# based on https://github.com/matplotlib/matplotlib/blob/v3.9.1/doc/conf.py#L757
repo_display_prefix = None # for viewcode instead
repo_display_prefix = f"{repo_url_prefix}/-/tree/{repo_tag}/testing/python_package"
if repo_display_prefix:
	import inspect

	extensions.append('sphinx.ext.linkcode')
	
	def linkcode_resolve(domain, info):
		"""
		Determine the URL corresponding to Python object
		"""
		if domain != 'py':
			raise ValueError(f'domain {domain}')
			# return None

		modname = info['module']
		fullname = info['fullname']
		# print(f'mof {modname} {fullname}')

		submod = sys.modules.get(modname)
		if submod is None:
			return None

		obj = submod
		for part in fullname.split('.'):
			try:
				obj = getattr(obj, part)
			except AttributeError:
				return None

		if inspect.isfunction(obj):
			obj = inspect.unwrap(obj)
		try:
			fn = inspect.getsourcefile(obj)
		except TypeError: 
			fn = None
		if not fn or fn.endswith('__init__.py'):
			try:
				fn = inspect.getsourcefile(sys.modules[obj.__module__])
			except (TypeError, AttributeError, KeyError):
				fn = None
		
		if not fn:
			return None

		try:
			source, lineno = inspect.getsourcelines(obj)
		except (OSError, TypeError):
			print(f'skipping{obj}...')
			lineno = None

		linespec = (f"#L{lineno:d}-L{lineno + len(source) - 1:d}"
					if lineno else "")

		startdir = Path(eden_simulator.__file__).parent.parent
		try:
			fn = os.path.relpath(fn, start=startdir).replace(os.path.sep, '/')
		except ValueError:
			return None

		# if not fn.startswith(('eden_simulator/')):
		# 	return None
		
		return (f"{repo_display_prefix}/{fn}{linespec}")
else:
	extensions.append('sphinx.ext.viewcode')


# nbsphinx

# hardcoded thumbnails for notebooks which normally don't show them, and for .gif thumbnails
nbsphinx_thumbnails = {
	'intro_neuroml': '_static/thumb_intro_neuroml.png',
	'intro_spatial': '_static/thumb_intro_spatial.png',
	'intro_lems'   : '_static/thumb_intro_lems.png'   ,
	'tut_net': '_static/thumb_tut_net.gif',
	'exa_lfp': '_static/thumb_exa_lfp.png',
	'extension_customsetup': '_static/thumb_extension_customsetup.png',
	'example_spatial_customsetup': '_static/thumb_example_spatial_customsetup.png',
	'extension_pointers'   : '_static/thumb_extension_pointers.png'   ,
	'extension_io'         : '_static/thumb_extension_io.png'         ,
	'extension_writable'   : '_static/thumb_extension_writable.png'   ,
	'extension_multiflux'  : '_static/thumb_extension_multiflux.png'  ,
	'example_pong': '_static/thumb_example_pong.png',
}

# This is processed by Jinja2 and inserted before each notebook https://github.com/spatialaudio/nbsphinx/blob/0.9.3/doc/conf.py#L43
# https://github.com/spatialaudio/nbsphinx/issues/419#issuecomment-603522609
# NEXT labpath instead of filepath?
# NEXT https://www.sphinx-doc.org/en/master/usage/theming.html#builtin-themes
# navigation_with_keys  body_max_width https://pydata-sphinx-theme.readthedocs.io/en/stable/user_guide/layout.html#references

# Add links to Binder, deep, colab
nbsphinx_prolog = r"""
{% set docname = 'docs/' + env.doc2path(env.docname, base=None) %}
{% set v = env.config.html_context %}

.. raw:: html

    <div id="run-notebook-online" class="admonition tip">
      This page was generated from
      <a class="reference external" href="{{ v.repo_url_prefix|e }}/blob/{{ v.repo_tag|e }}/{{ docname|e }}">{{ docname|e }}</a>.
      <span style="white-space: nowrap;">Interactive online version:</span>
	  
	  <span style="white-space: nowrap;"><a href="https://colab.research.google.com/github/{{ v.repo_github_user }}/{{ v.repo_name }}/blob/{{ v.repo_tag|e }}/{{ docname|e }}"><img alt="Colab badge" src="https://colab.research.google.com/assets/colab-badge.svg" style="vertical-align:text-bottom"></a></span>
	  
	  <span style="white-space: nowrap;"><a href="https://deepnote.com/launch?url={{ v.repo_url_prefix_github|e }}/blob/{{ v.repo_tag|e }}/{{ docname|e }}"><img alt="Deepnote badge" src="https://img.shields.io/badge/launch-deepnote-3793EF?logo=Deepnote" style="vertical-align:text-bottom"></a>.</span>
	  
	  <span style="white-space: nowrap;"><a href="https://mybinder.org/v2/gl/{{ v.repo_gitlab_user_binder|e }}%2F{{ v.repo_name }}/{{ v.repo_tag|e }}?filepath={{ docname|e }}"><img alt="Binder badge" src="https://mybinder.org/badge_logo.svg" style="vertical-align:text-bottom"></a></span>
	  
      <span style="white-space: nowrap;"><a href="{{ env.docname.split('/')|last|e + '.ipynb' }}" class="reference download internal" download>Download notebook</a>.</span>
    </div>

.. 
	raw:: latex
    \nbsphinxstartnotebook{\scriptsize\noindent\strut
    \textcolor{gray}{The following section was generated from
    \sphinxcode{\sphinxupquote{\strut {{ docname | escape_latex }}}} \dotfill}}
"""

# This is processed by Jinja2 and inserted after each notebook
nbsphinx_epilog = r"""
{% set docname = 'docs/' + env.doc2path(env.docname, base=None) %}
.. 
	raw:: latex
	%Perhaps clearpage or something? LATER


    \nbsphinxstopnotebook{\scriptsize\noindent\strut
    \textcolor{gray}{\dotfill\ \sphinxcode{\sphinxupquote{\strut
    {{ docname | escape_latex }}}} ends here.}}
"""
# Suppress prompt counter for publishing
# https://nbsphinx.readthedocs.io/en/0.9.5/configuration.html#nbsphinx_input_prompt
nbsphinx_input_prompt = '%.0s' 
nbsphinx_output_prompt = '%.0s'

# Allow funny content in %writefile cells
# see also https://github.com/spatialaudio/nbsphinx/issues/670, tried to override language_info.pygments_lexer but it didn't work
suppress_warnings += ['misc.highlighting_failure']

# common

# nitpicky = True LATER, complains about xyz, optional in 'returns' of docstrings

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
		get_more_assets(app)
	
	app.connect('builder-inited', builder_inited)
	
	def build_finished(app, exception):
		if exception is None:
			srcdir, outdir =  app.builder.srcdir, app.builder.outdir
			if app.builder.format == 'latex':
				# Grab the generated images and put them in the output's top level
				import shutil
				shutil.copytree(srcdir+'/_static', outdir+'/_static', dirs_exist_ok=True)
				# And delete this annoying log file for some reason, ir overwriting the out dir
				# shutil.rmtree(outdir+'/'+'eden-simulator'+'.ilg', ignore_errors=True)
			
	app.connect('build-finished', build_finished)
		

# -- Options for EPUB output --

# These are just defined to avoid Sphinx warnings related to EPUB: (from nbsphinx repo)
version = release
show_warning_types = True
suppress_warnings += ['epub.unknown_project_files']

# -- Options for HTML output --
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_static_path = ['_static']
html_theme = 'alabaster'
# html_css_files = ['_static/custom.css'] # HOTE: this does not seem to work, the css file is not even copied. But at least custom.css works for alabaster and pydata theme

# html_logo = "_static/eden_logo_white_bg.png"
# html_favicon = '_static/favicon.png' # NEXT https://sphinx-favicon.readthedocs.io/en/latest/quickstart.html#quickstart
html_theme_options = {
	'logo': 'eden_logo_white_bg.png', # _static is implicit?
}


# NEXT html_sidebars and nml link ! https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-html_sidebars

# Don't add .txt suffix to source files:
html_sourcelink_suffix = ''

# -- Options for LaTeX output --

# See https://www.sphinx-doc.org/en/master/latex.html
# more: https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-latex-output
# NEXT more: https://github.com/spatialaudio/nbsphinx/blob/0.9.4/doc/conf.py#L141
# NEXT display with env dependent behaviour...
# NEXT run linkcheck, can it also run on readthedocs https://github.com/spatialaudio/nbsphinx/blob/0.9.4/doc/conf.py#L115

# gotta love pygments, retrieve the macros because sphinx doesn't have the same ones apparently
import pygments; from pygments.formatters.latex import STYLE_TEMPLATE
pygments_latex_defs = STYLE_TEMPLATE % {'cp': 'PY', 'styles': ''}
# still doesn't look good, oh well, send a report to nbsphinx LATER

latex_elements = {
    'papersize': 'b5paper',
    # 'printindex': '',
#     'sphinxsetup': r"""
# HeaderFamily=\rmfamily\bfseries,
# div.note_border-TeXcolor={HTML}{E0E0E0},
# div.note_border-width=0.5pt,
# div.note_box-decoration-break=slice,
# div.warning_border-TeXcolor={HTML}{E0E0E0},
# div.warning_border-width=1.5pt,
# div.warning_background-TeXcolor={HTML}{FBFBFB},
# div.warning_box-decoration-break=slice,
# div.topic_box-shadow=none,
# div.topic_border-TeXcolor={HTML}{E0E0E0},
# div.topic_border-width=0.5pt,
# div.topic_box-decoration-break=slice,
# """,

#     'fontpkg': r"""
# \usepackage{mathpazo}
# \linespread{1.05}  % see http://www.tug.dk/FontCatalogue/urwpalladio/
# \setmainfont{TeX Gyre Pagella}[Numbers=OldStyle]
# \setmonofont{Latin Modern Mono Light}[Numbers=Lining]
# """,
    'preamble': r"""
\urlstyle{tt}

\usepackage{mathtools}
\usepackage[export]{adjustbox}

% for font fallback
\usepackage{fontspec, newunicodechar}
\newfontfamily{\fallbackfont}{Latin Modern Math}[Scale=MatchLowercase]
\DeclareTextFontCommand{\textfallback}{\fallbackfont}
% \newunicodechar{𝚛}{\textfallback{𝚛}}

\newcommand{\fallupchar}[2][\textfallback]{%
    \textup{#2}{#1{#2}}%
}

\newunicodechar{❗}{!} % LATER find a font for emojis
\newunicodechar{⠀}{ }

% \fallupchar{𝚛} LATER? missing \begin document?

\renewcommand{\sphinxsamedocref}[1]{#1}
\renewcommand{\sphinxcrossref}[1]{#1}
\renewcommand{\sphinxtermref}[1]{#1}

""" # + pygments_latex_defs
}
'''

\setmainfont{PT Serif}
\newfontfamily{\fallbackfont}{Linux Libertine O}[Scale=MatchLowercase]
\DeclareTextFontCommand{\textfallback}{\fallbackfont}
\newunicodechar{ɔ}{\textfallback{ɔ}}
\newunicodechar{ϱ}{\textfallback{ϱ}}

\begin{document}
Hellɔ woϱld.
\end{document}
'''
latex_toplevel_sectioning = 'part'

# LATER consider svg figures in pdf and even everywhere
# 	what about gouraud https://nbsphinx.readthedocs.io/en/0.9.5/code-cells.html#Plots
# 	and conf.py backends for mpl? or simply crank up the dpi?
# LATER more organised citations https://nbsphinx.readthedocs.io/en/0.9.4/a-normal-rst-file.html#citations
# though it would be best to convert doi url's into them, so that standalone notebooks can point to the doi url's
# -> for citations: just fetch all the hyperlinks and remap them to bib referenced with the help of a dict from url (or doi url) to bibtex item, it looks like the most sane *and* robust option.  Should work for html, ipynb and latex seamlessly.

latex_engine = 'lualatex'
latex_use_xindy = False

latex_show_urls = 'footnote'
latex_show_pagerefs = True # needed for print but just annoying for ebook ! LATER

latex_logo = '_static/eden_logo_white_bg.png'
# LATER perhaps reduce the fontsize of a bit, and convert unicore monospace to regular monospace.

linkcheck_ignore = [
	# JS based anchors
	'https://eden-simulator.org/repo#',
	'https://gitlab.com/c7859/neurocomputing-lab/Inferior_OliveEMC/eden/#',
	r'https://github.com\.*#L\.*',
	r'https://gitlab.com\.*#L\.*',
	# Flaky websites
	'https://v1.opensourcebrain.org/projects',
	'http://neuroml-db.org','https://neuroml-db.org',
	'https://onlinelibrary.wiley.com',
	
]
linkcheck_allowed_redirects = {
    # All HTTP redirections from the source URI to
    # the canonical URI will be treated as "working".
    r'https://sphinx-doc\.org/.*': r'https://sphinx-doc\.org/en/master/.*',
    r'https://doi.org/.*': r'.*',
    r'https://brian2\.readthedocs\.io/.*': r'https://brian2\.readthedocs\.io/en/stable/.*',
    r'https://nrn\.readthedocs\.io/.*': r'https://nrn\.readthedocs\.io/en/.*',
    r'https://docs\.pyvista\.org/api/core/_autosummary/pyvista\..*': r'https://docs\.pyvista\.org/api/core/_autosummary/pyvista\..*',
    r'https://gitlab\.com': r'https://about\.gitlab\.com',
}

# NEXT organise the structure better, with less inline toc and more intro pages...
# structure should be: intro (+ about i think!), user's guide, hacker's guide, python ref either at the end of the whole or at the end of the user's guide.

# done!
