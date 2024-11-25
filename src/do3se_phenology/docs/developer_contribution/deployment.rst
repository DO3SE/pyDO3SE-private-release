==================================
Deployment
==================================

The deployment process is still a WIP!

There are multiple deployment outputs to this project that need to be taken into account:

 - Compiled Python library
 - Compiled Python CLI
 - Documentation (Sphinx)
 - Python GUI


All of these should be automatically built using continuous deployment using CircleCI.
The circle CI config file is located in `.cicleci/config.yml`

All deployment should be proceeded by a test run to ensure the code still passes all tests
