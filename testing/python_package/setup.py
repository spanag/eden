import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="eden_tools",
    version="0.0.0",
    author="Sotirios Panagotou",
    author_email="sotiriospanagiotou0@gmail.com",
    description="Basic Python wrapper for EDEN simulator, and verification tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Environment :: Console",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Development Status :: 3 - Alpha",
        'Programming Language :: C++',
        "Programming Language :: Python :: 3",
    ],
    install_requires = [
        'pyNeuroML',
    ],
    python_requires='>=3.2',
)
