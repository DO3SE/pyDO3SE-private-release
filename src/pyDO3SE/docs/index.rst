.. pyDO3SE documentation master file, created by
   sphinx-quickstart on Wed Aug  5 10:57:16 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyDO3SE's documentation!
===================================


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   general/index
   end_user/index
   developer_contribution/index
   developer_external/index
   academic/index
   testing/do3se_doc1/index
   source/modules

**BROKEN**
Download a pdf copy of this documentation (WIP) :download:`pdf <pyDO3SE.pdf>`


Sections
--------

General Guidance
^^^^^^^^^^^^^^^^
This contains information that does not fit into the other categories such as legacy info.

End User Documentation
^^^^^^^^^^^^^^^^^^^^^^
Documentation on using the pyDO3SE Graphical User Interface (GUI) and Command Line
Interface (CLI)

Developer Contribution Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This contains guidance for developers working directly on the pyDO3SE code

External Developer Guidance
^^^^^^^^^^^^^^^^^^^^^^^^^^^
This contains guidance for developers that intend to use part or all of the pyDO3SE
model in their own model.

External Developer Guidance
^^^^^^^^^^^^^^^^^^^^^^^^^^^
This contains guidance for developers that intend to use part or all of the pyDO3SE
model in their own model.


pyDO3SE
^^^^^^^
This is the compiled source code documentation taken from docstrings in the source code.
This should be the ground truth for what is implemented in the model.

Download the executable
^^^^^^^^^^^^^^^^^^^^^^^
.. raw:: html

   <a href="./pyDO3SE.exe">pyDO3SE.exe</a>


Contributing to this documentation:
-----------------------------------
To contribute to this documentation first determine which section you need
to contribute to.
The documentation uses .rst files which are then built using Sphinx.
We also autoconvert word files(See below) and autogenerate a pdf (See below).
When changes have been made and pushed to the git repository these are autobuilt by CircleCi.
CircleCi then pushes the built documentation to the gh-pages branch.
The built documentation can then be viewed at https://sei-do3se.github.io/pyDO3SE/



Note: We have modified the sphinx default makefile to add production building.


.. include:: common/using_docx_files.rst
