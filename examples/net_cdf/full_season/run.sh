PROJECT_DIR="examples/net_cdf/full_season"
RUN_ID="demo_cli_run_merged"
python -m pyDO3SE run grid-run $PROJECT_DIR \
--multi-file-netcdf \
--start-input-index=6 \
--init-model \
--runid $RUN_ID \
--output-fields "td,td_dd,td_v,npp,dd,hr,dvi,f_phen,leaf_f_phen,canopy_height,lai,A_n,gsto,total_emerged_leaves" \
--netcdf-concat-dim "Time" \
--netcdf-regex-multi-file-filter '[0-9]{4}'


# PROJECT_DIR="examples/net_cdf/full_season"
# RUN_ID="demo_cli_run"
# python -m pyDO3SE run grid-run $PROJECT_DIR \
# --init-model \
# --runid $RUN_ID \
# --output-fields "pody,lai" \
# --netcdf-concat-dim "Time"