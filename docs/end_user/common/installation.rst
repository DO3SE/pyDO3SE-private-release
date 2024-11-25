=========================
Installation Instructions
=========================

Documentation on using the pyDO3SE model Graphical User Interface (GUI) and Command Line Interface (CLI).
Note: Most users will currently use the CLI.


Method 1 - Using the .exe file
==============================

This method involves using the .exe file located here:

.. raw:: html

   <a href="/pyDO3SE.exe">pyDO3SE.exe</a>

Note this will only work on Windows.
To run the exe file open windows powershell and navigate to the location of the pyDO3SE.exe file and run `./pyDO3SE.exe --help`.


Method 2 - Installing pyDO3SE into a Python environment
=======================================================

This method is best if you only need to run the model.

Setup ssh with github using the instructions here: https://github.com/sbland/python_for_modellers_guides/blob/main/using_version_control.md#ssh
Follow the instructions here on setting up a Python environment: https://docs.python.org/3/library/venv.html
Activate the environment then run:

:code:`pip install git+ssh://git@github.com/SEI-DO3SE/pyDO3SE@RELEASE`


Method 3 - Pulling the Git repository
=====================================

This method is best if you may need to make modifications to the code.
Setup version control using the instructions here: https://github.com/sbland/python_for_modellers_guides/blob/main/using_version_control.md#sphinx-quickstart

Clone the pyDO3SE repo with :code:`git clone git@github.com:SEI-DO3SE/pyDO3SE.git`
