==========
Single Run
==========

This page contains info on running single file DO3SE Cli runs.

You should use this run type if you want to run a single config file on a single
input data file.

Model Setup
===========

The DO3SE model requires the following input files:
 - A config file (config.json) :ref:`DO3SE Config Setup`
 - An input data file (input.csv) :ref:`Model Input Data Requirements`

You can include the following optional files
 - A base config file (base_config.json)

Running pyDO3SE with Single Run Setup
=====================================

For more info on running a single config run enter the following command into the
terminal: `python -m pyDO3SE run single --help`

Example Single Run
------------------
.. code-block:: bash
  :linenos:

  PROJECT_DIR=myprojectdir
  CONFIG_FILE=$PROJECT_DIR/configs/bangor_wheat.json
  DATA_FILE=$PROJECT_DIR/inputs/bangor_2015_hb_ww.csv
  OUTPUT_DIR=$PROJECT_DIR/runs/1/
  BASE_CONFIG_FILE=$PROJECT_DIR/base_config.json

  python -m pyDO3SE run single $CONFIG_FILE $DATA_FILE $OUTPUT_DIR -v \
  --base_config_file $BASE_CONFIG_FILE \
  --plot-fields "canopy_lai,canopy_height,dvi"


For more info see the :func:`.single` function.
