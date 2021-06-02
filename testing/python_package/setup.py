import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="eden_tools",
    version="0.0.0",
    author="Sotirios Panagoptou",
    author_email="spanagiotou@mail.ntua.gr",
    description="Basic Pyhton package for EDEN simulator, and testing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.2',
)
