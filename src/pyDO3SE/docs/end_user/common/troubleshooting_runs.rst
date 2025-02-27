====================
Troubleshooting Runs
====================

Documentation for troubleshooting model runs. This includes information on how to
understand the outputs, logs and to obtain more information on the model run.


Saving model in progress model state
====================================
You can save in progress model state at set intervals using the following CLI parameter
`dump_state_n=<interval>` where `n` is the row interval.
For example, to save the model state every 24 rows, you would use `dump_state_n=24`.


Example Cli

```bash
python -m pyDO3SE run batch \
$RUN_DIR \
--runid=$RUN_ID \
dump_state_n=960

```