.. EDEN documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

================================
Welcome to EDEN's documentation!
================================

.. include:: foreword.rst


.. raw:: latex
	
	% FOREWORD END
	% https://tex.stackexchange.com/questions/395856/switching-tocdepth-in-the-middle-of-a-document
	\newcommand{\changelocaltocdepth}[1]{%
	\addtocontents{toc}{\protect\setcounter{tocdepth}{#1}}%
	\setcounter{tocdepth}{#1}%
	}
	\newcommand{\changelocalsecnumdepth}[1]{%
	\addtocontents{toc}{\protect\setcounter{secnumdepth}{#1}}%
	\setcounter{secnumdepth}{#1}%
	}
	
	%\renewcommand{\partname}{Part}
	\renewcommand{\thepart}{\Alph{part}}

	
	\NewCommandCopy{\oldpart}{\part}
	\NewCommandCopy{\oldchapter}{\chapter}
	\renewcommand{\part}[1]{} % \textbf{\oldemph{#1}}
	
	\setcounter{secnumdepth}{-1}
	
	%\changelocaltocdepth{0}
	
	%\renewcommand{\chapter}[1]{\oldchapter*{#1}
	%\addcontentsline{toc}{chapter}{#1}
	%}

.. toctree::
	
	intro_latex

.. raw:: latex
	
	\setcounter{chapter}{0}
	\changelocalsecnumdepth{2}
	\changelocaltocdepth{2}
	
	\renewcommand{\part}{\oldpart}
	\renewcommand{\chapter}{\oldchapter}

.. toctree::

	user_guide_toctree

.. raw:: latex
	
	\setcounter{secnumdepth}{-2}
	
	%\renewcommand{\part}[1]{\oldpart*{#1}} 

.. toctree::
	
	refer_latex


