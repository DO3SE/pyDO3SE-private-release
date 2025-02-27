CONFIG_FILE=examples/bangor_2015_multilayer/configs/bangor_wheat.json
DATA_FILE=examples/bangor_2015_multilayer/inputs/bangor_2015_hb_ww.csv
OUTPUT_DIR=examples/bangor_2015_multilayer/runs/from_cli_d/
BASE_CONFIG_FILE=examples/bangor_2015_multilayer/base_config.json

python -m pyDO3SE run single $CONFIG_FILE $DATA_FILE $OUTPUT_DIR "debug=True" \
--base_config_file $BASE_CONFIG_FILE \
--plot-fields "fO3_l,f_LS,f_LA,fst_acc"

# python -m pyDO3SE analysis compare-outputs examples/bangor_2015_multilayer/runs/from_cli_b/ tmp/comparisons --additional-input-dirs examples/bangor_2015_multilayer/runs/from_cli/