====================
pyDO3SE CLI Guidance
====================

Guidance on using the pyDO3SE CLI


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   grid_run


Using the CLI
=============
To run the pyDO3SE CLI
- Follow the instructions for installing pyDO3SE.
- Run `python -m pyDO3SE --help`
- Follow the instructions shown to run use the cli.



.. click:: pyDO3SE.tools.cli.root:cli
   :nested: full

.. click:: pyDO3SE.tools.cli.root:cli
   :prog: pyDO3SE
   :nested: full


Performing a multi config run
=============================

To perform a multi config run `python -m pyDO3SE run batch --help`