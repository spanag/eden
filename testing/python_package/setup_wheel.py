

# Before loading setuptools, play with cmdline args here. Because setuptools will complain for the parameters it doesn't recognize
import sys
import argparse

# Add custom parameters to setup.py
parser = argparse.ArgumentParser( allow_abbrev=False )
parser.add_argument(
    '--package-version',
    help='Version tag of the package to be built',
    required=True
)
# boolean optional argparse parameters: https://stackoverflow.com/a/63085690
parser.add_argument(
    '--no-eden-exe',
    help='Exclude binary executables from package',
    action="store_true" # false by default
)

(config_options, argv) = parser.parse_known_args()
sys.argv = [sys.argv[0]] + argv  # Write the remaining arguments back to `sys.argv` for distutils to read
assert (config_options.package_version)

# TODO maybe write version.py file here...
package_name = "eden_simulator"
package_dir = package_name


exe_files_location = 'data/bin'

# Now, setuptools can be loaded safely
import setuptools
import wheel
import distutils
import sysconfig


# load README
with open("README_wheel.md", "r") as fh:
    long_description = fh.read()

# Restrict tag to platform dependent, without restricting abi, and don't actually build any extension
# As seen on https://github.com/Yelp/dumb-init/blob/master/setup.py#L11-L27
# try:
from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
class platform_specific_bdist_wheel(_bdist_wheel):

    def finalize_options(self):
        _bdist_wheel.finalize_options(self)
        # Mark us as not a pure python package
        self.root_is_pure = False
    
    # Tag elements are: python version, abi, platform, ...
    # No abi requirement, keep the tags
    def get_tag(self):
        return ('py3', 'none',) + _bdist_wheel.get_tag(self)[2:]
        
# except None: # ImportError:
#     platform_specific_bdist_wheel = None # well the wheel package is really necessary...

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

def GetBuildScriptsExcept(list_of_binaries):
    
    import distutils.command.build_scripts
    
    class my_build_scripts(distutils.command.build_scripts.build_scripts):
        """Copy each script listed in 'self.scripts' *as is*, without
        attempting to adjust the !# shebang. The default build_scripts command
        in python3 calls tokenize to detect the text encoding, treating all
        scripts as python scripts. But all our 'scripts' are native C
        executables, thus the python3 tokenize module fails with SyntaxError
        on them. Here we just skip the if branch where distutils attempts
        to adjust the shebang.
        
        Modified to apply only to pre-specified binaries, otherwise the typical processing applies.
        """
        
        def copy_scripts(self):
            all_scripts = self.scripts
            
            binary_scripts = [ x for x in all_scripts if x     in list_of_binaries]
            othery_scripts = [ x for x in all_scripts if x not in list_of_binaries]
            
            self.scripts = binary_scripts
            (b_outfiles, b_updated_files) = self.copy_binary_scripts()
            
            self.scripts = othery_scripts
            (o_outfiles, o_updated_files) = super().copy_scripts()
            
            self.scripts = all_scripts
            return (b_outfiles + o_outfiles, b_updated_files+o_updated_files )
        
        def copy_binary_scripts(self):
            self.mkpath(self.build_dir)
            outfiles = []
            updated_files = []
            for script in self.scripts:
                # print('*** '+ script)
                self._copy_script_binary(script, outfiles, updated_files)
                
            return outfiles, updated_files
        '''
        def _copy_script(self, script, outfiles, updated_files):
            if script not in list_of_binaries:
                super()._copy_script(script, outfiles, updated_files)
            else:
                self._copy_script_binary(script, outfiles, updated_files)
        '''
        def _copy_script_binary(self, script, outfiles, updated_files):
            import os
            from distutils.util import convert_path
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
                    return

            if f:
                f.close()
            
            updated_files.append(outfile)
            self.copy_file(script, outfile)
            
    return my_build_scripts

# Executable "console scripts" are either Unix-style #! executable scripts, or binary files. 
# Windows binaries are greated through an ingenious hack that appends the Python command to a 'launcher' binary, which then finds what's stuck on top of it and runs it properly. See https://github.com/pypa/distlib/blob/0.3.6/distlib/scripts.py#L262 and https://github.com/vsajip/simple_launcher/blob/master/launcher.c
# TODO Therefore we can probably eliminate messing with the exe "scripts" in our case, by using console scripts only that wrap around the executable(s) placed on data/ anyway.
# The "override stuff in distutils" technique still applies when the "script" really must be the application binary.

def _get_scripts():
    import platform
    
    python_scripts = []
    binary_scripts = []
    if not config_options.no_eden_exe: 
        exe_script_names = [
            'eden',
        ]
        if platform.system() == 'Windows':
            extension = '.exe'
        else:
            extension = ''
            
        # if sys.platform == "darwin":
            # return [] # TODO a python wrapper instead
            
        binary_scripts = [package_name+'/'+exe_files_location+'/%s%s' % (script_name, extension)
            for script_name in exe_script_names]
        
    return python_scripts, binary_scripts

[python_scripts, binary_scripts] = _get_scripts()
scripts_list = python_scripts + binary_scripts

cmd_my_build_scripts = GetBuildScriptsExcept(binary_scripts)

if config_options.no_eden_exe:
    # Python only wheel, in case eden is installed and in PATH by other means
    cmd_my_bdist_wheel = wheel.bdist_wheel.bdist_wheel
    cmd_my_install = distutils.command.install.install
    my_distclass = setuptools.dist.Distribution
    zip_safe = True
    package_data_locations = []
    
    # cmd_my_bdist_wheel =  platform_specific_bdist_wheel
    # cmd_my_install = InstallPlatlib
    # my_distclass = BinaryDistribution
    zip_safe = False
else:
    # Wheel contains compiled code and is thus platform specific
    cmd_my_bdist_wheel =  platform_specific_bdist_wheel
    cmd_my_install = InstallPlatlib
    my_distclass = BinaryDistribution
    zip_safe = False
    package_data_locations = [exe_files_location+'/*']




setuptools.setup(
    name=package_name,
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
    packages=[package_name],
    # for future reference: whether binaries are included or silently dropped is up to the python packaging system's whims, the whims of the day were last tracked here https://github.com/pypa/setuptools/issues/3340#issuecomment-1219383976
    # packages=setuptools.find_namespace_packages(where=package_dir),
    # packages=setuptools.find_namespace_packages(where="."),
    # package_dir={"":package_dir},
    # package_dir={"":'.'},
    # package_dir={package_dir:package_dir},
    package_data={package_name: package_data_locations},
    include_package_data=False,
    # scripts 
    scripts=scripts_list,  
    # scripts=['eden'], 
    install_requires = [
        'pyneuroml',
        # 'setuptools', # due to customised setup step ... but it should already be in place to install the wheel right?
        'lxml',
    ] + (['h5py <= 2.10.*'] if (sysconfig.get_platform() == 'win32') else []) # h5py wheels are missing since, and pip doesn't know that h5py source is tough to build
    ,

    python_requires='>=3.2',
    platforms=['windows','mac','linux'],
    cmdclass={
        'bdist_wheel': cmd_my_bdist_wheel,
        'build_scripts': cmd_my_build_scripts,
        'install': cmd_my_install,
    },
    distclass=my_distclass,
    zip_safe=zip_safe,
)
