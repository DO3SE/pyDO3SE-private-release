.. do3se_phenology documentation master file, created by
   sphinx-quickstart on Wed Aug  5 10:57:16 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to do3se_phenology's documentation!
===========================================


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   end_user/index
   developer_contribution/index
   source/modules

**BROKEN**
Download a pdf copy of this documentation (WIP) :download:`pdf <do3se_phenology.pdf>`


Sections
--------

_End User Documentation
^^^^^^^^^^^^^^^^^^^^^^^
Documentation on using the do3se_phenology Graphical User Interface (GUI) and Command Line
Interface (CLI)

_Developer Contribution Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This contains guidance for developers working directly on the do3se_phenology code


_External Developer Guidance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This contains guidance for developers that intend to use part or all of the do3se_phenology
model in their own model.


_do3se_phenology Source
^^^^^^^^^^^^^^^^^^^^^^^
This is the compiled source code documentation taken from docstrings in the source code.
This should be the ground truth for what is implemented in the model.


Contributing to this documentation:
-----------------------------------
To contribute to this documentation first determine which section you need
to contribute to.
The documentation uses .rst files which are then built using Sphinx.
We also autoconvert word files(See below) and autogenerate a pdf (See below).
When changes have been made and pushed to the git repository these are autobuilt by CircleCi.
CircleCi then pushes the built documentation to the gh-pages branch.




Note: We have modified the sphinx default makefile to add production building.
