# The DO3SE Crop Impact Assessment Tool

DO3SE (Deposition of Ozone for Stomatal Exchange) is a dry deposition model designed to estimate the total and stomatal deposition (or flux) of ozone (O3) to selected European land-cover types and plant species.

# About pyDO3SE2

pyDO3SE2 is a full rebuild of the original pyDO3SE model that was written in Fortran.

# Full documentation

Full documentation is available at

https://pydo3se-docs.onrender.com/
This is compiled from the docs directory using sphinx and built on circle ci

# Using pyDO3SE End User guide

## Environment Setup

1. Setup Python 3.8 by following the instructions here: https://github.com/sbland/python_for_modellers_guides/blob/main/setup_python_environment.md
2. Make sure to create and activate a virtual environment.
3. Install _wheel_ to help build the packages with `pip install wheel`.
4. Install pyDO3SE with `pip install git+ssh://git@github.com/SEI-DO3SE/pyDO3SE@RELEASE`
5. Access the CLI with `python -m pyDO3SE --help`

To update pyDO3SE run `pip install git+ssh://git@github.com/SEI-DO3SE/pyDO3SE@RELEASE --upgrade`
Use `pip install git+ssh://git@github.com/SEI-DO3SE/pyDO3SE@RELEASE --upgrade` to quickly update pyDO3SE only.

## Running the model

### CLI

To run the model from the cli (cmd or bash)

1. activate the virtual environment
2. run `python -m pyDO3SE.py` and follow the instructions

To run a file run `python -m pyDO3SE.py run single --help` for instructions.

### Config Instructions

To get an example config file run `python -m pyDO3SE.py config generate --help` and follow the instructions.
Additional guidance on the requirements of the config file is in progress. In the mean time refer to the files
`Config_Shape.py` and `ConfigLandCover.py`.

### GUI (NOT BEING DEVELOPED CURRENTLY)

The gui is currently a WIP and uses `streamlit` to run a simple web browser based ui.
For more info go to the README.md in the pyDO3SE_gui directory

To run the gui from the cli:

1. first install the `requirements/gui.txt` file
2. Run the gui with `streamlit run pyDO3SE_gui.py` the follow the instructions

---

# Config Files

Examples of the config files that define the inital parameters are located here:
`examples/spanish_wheat/spanish_wheat_config.json`

The file should be a json format which matches the model config located in `pyDO3SE/Config/Config_Shape.py`

# Multirun

To run multiple configs on multiple files use `python -m pyDO3SE.py run batch --help`

The files for a multirun should be located in a config and input directories. Each config will be ran for each
input file then saved into a output directory. The output directory will contain a directory for each config.

# Output comparison graphs

To quickly create comparison graphs of outputs run `python -m pyDO3SE.py analysis compare-outputs --help`.
This will need an input directory that contains output files; an output directory to save
the graphs too and then a list of fields to graph. Fields requested and the column headings in the output
must match those in `pyDO3SE/Output/Output_Shape.py` Use the first value in the field object.

# Contributing

## Updating version

The `bump2version` package is used for version management.
To update the version type `bumpversion patch|minor|major`.
This will update the version in python.py and version.py

## Environment Setup

Clone this repository by following the instructions here: https://github.com/sbland/python_for_modellers_guides/blob/main/using_version_control.md
Note: You will need to setup a SSH key.

Setup Python 3.8 by following the instructions here: https://github.com/sbland/python_for_modellers_guides/blob/main/setup_python_environment.md

Make sure to create and activate a virtual environment.

1. Complete all the prerequisites above
2. Install _wheel_ to help build the packages with `pip install wheel`.
3. Install from the requirements file using `pip install -r requirements/<requirementsfile>.txt`
   - **common.txt** - required by all environments
   - **development.txt** - required for contributing to the code (optional)
   - **docs.txt** - required to contribute to the sphynx documentation (optional)
   - **jupyter.txt** - required for working with the notebooks (optional)
4. Test everything worked by running `python pyDO3SE_cli.py demo`

---

## Running the test suite

The tests use pytest https://docs.pytest.org

To run the test suite run the following from the root of
the repo (i.e. next to this README file):

1. activate the virtual environment
2. `python setup.py test`

You will be given a summary of the tests that passed/failed along with a code
coverage report outlining which lines in the code base have not been tested by
the current test suite.
