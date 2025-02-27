=====================
Updating Instructions
=====================

Documentation on using the pyDO3SE model Graphical User Interface (GUI) and Command Line Interface (CLI).
Note: Most users will currently use the CLI.


Method 1 - Using the .exe file
==============================

If you are using the .exe file please redownload the latest file below:

.. raw:: html

   <a href="/pyDO3SE.exe">pyDO3SE.exe</a>

Note this will only work on Windows.
To run the exe file open windows powershell and navigate to the location of the pyDO3SE.exe file and run `./pyDO3SE.exe --help`.


Method 2 - Installing pyDO3SE into a Python environment
=======================================================

If you installed pyDO3SE in a Python Environment you will need to run the following to update the model

:code:`pip install -U git+ssh://git@github.com/SEI-DO3SE/pyDO3SE@RELEASE`


Method 3 - Pulling the Git repository
=====================================

If you cloned the pyDO3SE git repository you will need to pull the changes as below:

Clone the pyDO3SE repo with :code:`git pull`

If you have made any changes you will need to either stash these or commit them.
