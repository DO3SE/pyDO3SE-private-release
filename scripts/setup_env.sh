rm -R venv
python -m venv venv
source venv/bin/activate

pip install -r requirements/common.txt
pip install -e src/pyDO3SE
pip install -e src/thermal_time
pip install -e src/do3se_phenology
pip install -e src/do3se_met