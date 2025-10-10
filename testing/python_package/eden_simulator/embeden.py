'''Routines to access the package's embedded binaries, if they exist.'''

# NB this is safe ONLY because the package is NOT zip safe!!
# NB the only way to make the package zip safe is by extracting ALL bundled dlls along with the exe! 
# Perhaps delvewheel, auditwheel, delocate & friends might help in the future; still more dependencies that the binary knows about will have to be unpacked in the expected locations.
# Until then this package is NOT zip safe.

import sys, platform
exe_extension = ".exe" if platform.system() == 'Windows' else ""
eden_bundled_exe = "data/bin/eden"+exe_extension

# returns None if not found
def get_exe_path():
	eden_bundled_exe_filename = None
	if sys.version_info >= (3,9):
		# new way https://importlib-resources.readthedocs.io/en/latest/migration.html#pkg-resources-resource-filename
		from importlib import resources as importlib_resources
		ref = importlib_resources.files(__name__) / eden_bundled_exe
		if ref.is_file():
			# NB I won't do this the "safe" way because there are YET more files that should be unpacked alongside the exe!
			# See also the start of this file. The only cure is for the package to be always unpacked as a tree.
			with importlib_resources.as_file(ref) as path: # path would disappear after this scope, were the package compressed.
				eden_bundled_exe_filename = str(path)
	else:
		# old way
		import pkg_resources
		if pkg_resources.resource_exists(__name__, eden_bundled_exe):
			eden_bundled_exe_filename = pkg_resources.resource_filename(__name__, eden_bundled_exe)
		
	return eden_bundled_exe_filename
