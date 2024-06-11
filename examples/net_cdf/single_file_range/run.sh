PROJECT_DIR="examples/net_cdf/single_file_range"
RUN_ID="demo_cli_run_merge"
python -m pyDO3SE run grid-run $PROJECT_DIR \
--multi-file-netcdf \
--init-model \
--runid $RUN_ID \
--output-fields "pody,lai" \
--netcdf-concat-dim "Time" \
--netcdf-regex-multi-file-filter '[0-9]{4}'


PROJECT_DIR="examples/net_cdf/single_file_range"
RUN_ID="demo_cli_run"
python -m pyDO3SE run grid-run $PROJECT_DIR \
--init-model \
--runid $RUN_ID \
--output-fields "pody,lai" \
--netcdf-concat-dim "Time"