

# Before loading setuptools, play with cmdline args here. Because setuptools will complain for the parameters it doesn't recognize
import sys
import argparse
import os
# Add custom parameters to setup.py
parser = argparse.ArgumentParser( allow_abbrev=False )
parser.add_argument(
    '--package-version',
    help='Version tag of the package to be built',
    required=True
)

(config_options, argv) = parser.parse_known_args()
sys.argv = [sys.argv[0]] + argv  # Write the remaining arguments back to `sys.argv` for distutils to read
assert (config_options.package_version)

# Now, setuptools can be loaded safely
import setuptools

# load README
with open("README_wheel.md", "r") as fh:
    long_description = fh.read()

# Restrict tag to platform dependent, without restricting abi, and don't actually build any extension
# As seen on https://github.com/Yelp/dumb-init/blob/master/setup.py#L11-L27
try:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel

    class my_bdist_wheel(_bdist_wheel):

        def finalize_options(self):
            _bdist_wheel.finalize_options(self)
            # Mark us as not a pure python package
            self.root_is_pure = False
        
        # Tag elements are: python version, abi, platform, ...
        # No abi requirement, keep the tags
        def get_tag(self):
            return ('py3', 'none',) + _bdist_wheel.get_tag(self)[2:]
        
except ImportError:
    my_bdist_wheel = None

from setuptools.dist import Distribution
class BinaryDistribution(Distribution):
    def is_pure(self):
        return False

    def has_ext_modules(self):
        return True

# Actually copy the binaries as "scripts" of machine code
# TODO maybe it is wiser (and more space efficient) to make a true Python wrapper, if access to the true binary is not necessary (think of env vars, and other things that must be transparent, also think of the delay in launching python?)
# As seen on https://github.com/adobe-type-tools/afdko/blob/develop/setup.py#L55 copr. Adobe
import distutils.command.install as _base_install
class InstallPlatlib(_base_install.install):
    """This is to force installing all the modules to the non-pure, platform-
    specific lib directory, even though we haven't defined any 'ext_modules'.
    The distutils 'install' command, in 'finalize_options' method, picks
    either 'install_platlib' or 'install_purelib' based on whether the
    'self.distribution.ext_modules' list is not empty.
    Without this hack, auditwheel would flag the wheel as invalid since
    it contains native executables inside the pure lib folder.
    TODO Remove this hack if/when in the future we install extension modules.
    """

    def finalize_options(self):
        _base_install.install.finalize_options(self)
        self.install_lib = self.install_platlib

import distutils.command.build_scripts
class my_build_scripts(distutils.command.build_scripts.build_scripts):

    def copy_scripts(self):
        """Copy each script listed in 'self.scripts' *as is*, without
        attempting to adjust the !# shebang. The default build_scripts command
        in python3 calls tokenize to detect the text encoding, treating all
        scripts as python scripts. But all our 'scripts' are native C
        executables, thus the python3 tokenize module fails with SyntaxError
        on them. Here we just skip the if branch where distutils attempts
        to adjust the shebang.
        """
        import os
        from distutils.util import convert_path
        self.mkpath(self.build_dir)
        outfiles = []
        updated_files = []
        for script in self.scripts:
            script = convert_path(script)
            outfile = os.path.join(self.build_dir, os.path.basename(script))
            outfiles.append(outfile)

            try:
                f = open(script, "rb")
            except OSError:
                if not self.dry_run:
                    raise
                f = None
            else:
                first_line = f.readline()
                if not first_line:
                    f.close()
                    self.warn("afdko: %s is an empty file (skipping)" % script)
                    continue

            if f:
                f.close()
            updated_files.append(outfile)
            self.copy_file(script, outfile)

        return outfiles, updated_files

def _get_scripts():
    import platform
    
    script_names = [
        'eden',
    ]
    if platform.system() == 'Windows':
        extension = '.exe'
    else:
        extension = ''
        
    if sys.platform == "darwin":
        return [] # TODO a python wrapper instead
        
    scripts = ['bin/%s%s' % (script_name, extension)
        for script_name in script_names]
        
    return scripts

package_data_locations = ['../bin/*']

zip_safe = True

if sys.platform == "darwin":
    package_data_locations += ['../bin/dylibs/*']
    zip_safe = False
    
    # alternatively, for dylib hierarchy (TODO)
    # or maybe graft bin/ instead?
    # for root, _, filenames in os.walk('bin/dylibs'):
    #         for fname in filenames:
    #             fullname = join(root, fname)
    #             scripts = scripts + [ for file_name in os.walk('bin/dylibs')]

setuptools.setup(
    name="eden_simulator",
    version=config_options.package_version,
    description="Standalone Python wheels for the EDEN neural simulator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Sotirios Panagotou",
    author_email="s.panagiotou@erasmusmc.nl",
    url='https://gitlab.com/neurocomputing-lab/Inferior_OliveEMC/eden',
    download_url='https://gitlab.com/neurocomputing-lab/Inferior_OliveEMC/eden/-/archive/main/eden-main.zip',
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Environment :: Console",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS",
        "Operating System :: POSIX :: Linux",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Development Status :: 3 - Alpha",
        'Programming Language :: C++',
        "Programming Language :: Python :: 3",
    ],
    keywords=['simulator','simulation','HPC','neuroscience','NeuroML'],
    # project_urls
    packages=['eden_simulator'],
    package_data={'eden_simulator': package_data_locations},
    # include_package_data=True,
    # scripts 
    scripts=_get_scripts(),  
    # scripts=['eden'], 
    install_requires = [
        'pyneuroml',
    ],
    python_requires='>=3.2',
    platforms=['linux'],
    cmdclass={
        'bdist_wheel': my_bdist_wheel,
        'build_scripts': my_build_scripts,
        'install': InstallPlatlib,
    },
    distclass=BinaryDistribution,
    zip_safe=zip_safe,
)
