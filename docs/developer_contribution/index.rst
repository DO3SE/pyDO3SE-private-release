=================================
Developer Contribution Guidelines
=================================

This contains guidance for developers working directly on the pyDO3SE code. Also ensure you familiarise yourself with the end user documentation.


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   developer_environment
   testing
   deployment
   git_workflows
   functional_approach_to_state_models
   model_data
   model_notes


Linting
-------
To ensure a consistent style this project uses flake8 for linting.

pre-commit
----------
This project is set up with pre-commit to ensure nothing is commited without first passing some linting rules.
Ensure these are setup by running `pre-commit install`

testing
-------
Snapshot testing is used to snapshot the output of some tests that contain a lot of output data.

Docstrings
----------
The docstrings follow the numpy format.

Git Workflow
------------

- Master branch contains most up to date code
- Stable releases are tagged with their release number
- Additional branches are created for features that have not yet been fully implemented
- Branches should not be used to store alternate versions of the model unless they are assumed to be static releases that will not be updated.
This avoides trying to keep track of many different versions. Instead you should attempt to use the process runner to create an alternative set of default processes.

Python features
---------------
A few python features that are used that could be considered new or not always used.

Dataclasses
^^^^^^^^^^^
https://docs.python.org/3/library/dataclasses.html

Typing
^^^^^^
https://docs.python.org/3/library/typing.html


Process Runner (Proflow)
^^^^^^^^^^^^^^^^^^^^^^^^
In order to control the flow of processes and data in the model and make it easier to debug where state is
being changed we use the process runner library.
Instructions on using the library can be found here:
https://github.com/sbland/proFlow

The processes that are run by default are located in `pyDO3SE/Defaults/default_processes.py`


Directory Structure
-------------------

::

   root
   |-- docs/ => Contains the documentation
   |-- examples/ => contains example data and configs
   |-- helpers/ => A selection of python helper functions
   |-- notebooks/ => Jupyter notebooks for testing parts of the model
   |-- pyDO3SE/ => The pyDO3SE package
   |   |-- Config/ => Model config
   |   |-- constants/ => Model constants
   |   |-- Defaults/ => Default processes to run in the model
   |   |-- External_State/ => External data
   |   |-- Model_state/ => The model state
   |   |-- Parameters/ => DEPRECIATED
   |   |-- plugins/ => pyDO3SE modules/plugins
   |-- pyDO3SE_gui/ => pyDO3SE gui package
   |-- requirements/ => Contains the requirement.txt familiarise
   |-- test/ => functional tests
   |-- vendor/ => contains unreleased python packages that are required by pyDO3SE