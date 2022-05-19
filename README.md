# clamdb
CLaM Database Management

For a detailed description, please see [the documentation website](https://fuzzy-couscous-01244cb5.pages.github.io/)
# Installation
```
# Install the 'githubinstall' package from CRAN, if you have not already done so.
install.packages("githubinstall")

# Install clamr using githubinstall (clamdb dependency)
githubinstall::githubinstall("clamr")

# Install clamdb using githubinstall
githubinstall::githubinstall("clamdb")
```

# Additional Installation
`clamdb` has functionality to import and build `.msp` libraries, however, this requires
some additional installation.  Specifically, some `python` code that comes bundled with
this package must also be installed.

All of the python code `clamdb` requires is located in the `inst/python` folder.

To install this code, you must first install the following python packages:
```
pandas
rdkit
```
You must also have configured a conda environment, which is necessary for using the R package
[reticulate](https://rstudio.github.io/reticulate/articles/python_packages.html).

Once your conda environment has been configured, you can install the python code here by
cloning this repository, navigating to the `inst/python` subfolder, and executing
`python setup_mschem.py sdist bdist_wheel`, e.g.

```
git clone <this-repo> <your-system-path>
cd <your-system-path>/clamdb/inst/python
python setup_mschem.py sdist bdist_wheel
```
Where `your-system-path` would be specific to where you would like to clone this repository.

