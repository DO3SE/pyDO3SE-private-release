PROJECT_DIR="."
RUN_ID="demo_cli_run_merged2"
python -m pyDO3SE run grid-run $PROJECT_DIR \
--init-model \
--runid $RUN_ID \
--loglevel 2 \
--output-fields "td,td_dd,td_v,npp,dd,hr,dvi,f_phen,leaf_f_phen,canopy_height,lai,A_n,gsto,total_emerged_leaves"
