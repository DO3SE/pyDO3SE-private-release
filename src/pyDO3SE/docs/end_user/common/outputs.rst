======================
DO3SE End User Outputs
======================

Documentation on the DO3SE outputs.


The output directory will be populated as below::

- --- runs
- ----- <config_1>
- ------- <input_1>
- --------- config.json
- --------- final_state.json
- --------- notes.log
- --------- processed_config.json
- --------- external_data.csv
- --------- <config>_<input>_out.csv
- --------- <field_1>.png
- --------- <field_2>.png
- --------- ... additional field plots



config.json
    A copy of the input config file
final_state.json
    The final state of the model
notes.log
    Any additional notes on the model run including model version;
    time taken to run etc.
processed_config.json
    The actual model config that has been ran including any default,
    config processing etc.
external_data.csv
    The processed external data that has been used for the model run.
    This includes any precalculated values.
<config>_<input>_out.csv
    The output of the model run
<field_1>.png
    Plots of fields selected
output_fields_info.csv
    Meta data for each heading in the output file