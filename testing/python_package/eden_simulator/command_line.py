#!/usr/bin/python3
# see also: https://github.com/scikit-build/cmake-python-distributions/blob/3.26.3/src/cmake/__init__.py
import os
import subprocess
import sys

from . import embeden
eden_bundled_exe_filename = embeden.get_exe_path()

def _DoCall(binary_filename):
	raise SystemExit( subprocess.call([binary_filename] + sys.argv[1:], close_fds=False))

def runEden():
	if not eden_bundled_exe_filename:
		# then command line should not be called with 'eden' because that could be this script again, abort instead
		raise NotImplementedError("EDEN binary is not bundled!")
	
	return _DoCall(eden_bundled_exe_filename)
