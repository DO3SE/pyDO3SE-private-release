# Initial clone
PYDO3SE_LOCAL_PATH="src/pyDO3SE"
THERMAL_TIME_LOCAL_PATH="src/thermal_time"
DO3SE_PHENOLOGY_LOCAL_PATH="src/do3se_phenology"
DO3SE_MET_LOCAL_PATH="src/do3se_met"

PYDO3SE_REMOTE="git@github.com:SEI-DO3SE/pyDO3SE"
THERMAL_TIME_REMOTE="git@github.com:SEI-DO3SE/thermal_time"
DO3SE_PHENOLOGY_REMOTE="git@github.com:SEI-DO3SE/do3se_phenology"
DO3SE_MET_REMOTE="git@github.com:SEI-DO3SE/do3se_met"


git subtree add --prefix $PYDO3SE_LOCAL_PATH $PYDO3SE_REMOTE RELEASE --squash
git subtree add --prefix $THERMAL_TIME_LOCAL_PATH $THERMAL_TIME_REMOTE RELEASE --squash
git subtree add --prefix $DO3SE_PHENOLOGY_LOCAL_PATH $DO3SE_PHENOLOGY_REMOTE RELEASE --squash
git subtree add --prefix $DO3SE_MET_LOCAL_PATH $DO3SE_MET_REMOTE RELEASE --squash
