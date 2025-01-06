PROJECT_DIR="examples/soil_moisture"
OUTPUT_FIELDS="gsto,canopy_lai"
python -m pyDO3SE run batch $PROJECT_DIR \
--runid=1 \
--run_comparisons \
--output_fields=$OUTPUT_FIELDS \
--compare_fields=$OUTPUT_FIELDS
