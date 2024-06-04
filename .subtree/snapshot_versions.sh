PYDO3SE_LOCAL_PATH="src/pyDO3SE"
THERMAL_TIME_LOCAL_PATH="src/thermal_time"
DO3SE_PHENOLOGY_LOCAL_PATH="src/do3se_phenology"
DO3SE_MET_LOCAL_PATH="src/do3se_met"

PYDO3SE_VERSION=`cat $PYDO3SE_LOCAL_PATH/VERSION.txt`
THERMAL_TIME_VERSION=`cat $THERMAL_TIME_LOCAL_PATH/VERSION.txt`
DO3SE_PHENOLOGY_VERSION=`cat $DO3SE_PHENOLOGY_LOCAL_PATH/VERSION.txt`
DO3SE_MET_VERSION=`cat $DO3SE_MET_LOCAL_PATH/VERSION.txt`


echo "pyDO3SE version: $PYDO3SE_VERSION"
echo "thermal_time version: $THERMAL_TIME_VERSION"
echo "do3se_phenology version: $DO3SE_PHENOLOGY_VERSION"
echo "do3se_met version: $DO3SE_MET_VERSION"



# Store version in file named SUB_VERSIONS.txt
echo "pyDO3SE version: $PYDO3SE_VERSION" > SUB_VERSIONS.txt
echo "thermal_time version: $THERMAL_TIME_VERSION" >> SUB_VERSIONS.txt
echo "do3se_phenology version: $DO3SE_PHENOLOGY_VERSION" >> SUB_VERSIONS.txt
echo "do3se_met version: $DO3SE_MET_VERSION" >> SUB_VERSIONS.txt
# == End User Documentation == #