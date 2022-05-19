from setuptools import setup, find_packages

VERSION = '1.0.0' 
DESCRIPTION = 'mschem python package for use with clamdb'
LONG_DESCRIPTION = 'Python utilities from the rdkit package wrapped in convenience functions for clamdb.'

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="mschem", 
        version=VERSION,
        author="Sean Hackett",
        author_email="<sean@calicolabs.com>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=['rdkit', 'pandas'],
        keywords=['mschem', 'open_clam', 'clamdb'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Education",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
        ]
)