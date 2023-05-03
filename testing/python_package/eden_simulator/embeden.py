# Routines to access the package's embedded binaries, if they exist

import pkg_resources
import platform
exe_extension = ".exe" if platform.system() == 'Windows' else ""
eden_bundled_exe = "data/bin/eden"+exe_extension

# returns None if not found
def get_exe_path():
	eden_bundled_exe_filename = None
	if pkg_resources.resource_exists(__name__, eden_bundled_exe):
		eden_bundled_exe_filename = pkg_resources.resource_filename(__name__, eden_bundled_exe)
	
	return eden_bundled_exe_filename
