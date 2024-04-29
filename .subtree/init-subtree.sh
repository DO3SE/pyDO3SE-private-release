# Initial clone
PYDO3SE_LOCAL_PATH="src/pyDO3SE"
THERMAL_TIME_LOCAL_PATH="src/thermal_time"
DO3SE_PHENOLOGY_LOCAL_PATH="src/do3se_phenology"
DO3SE_MET_LOCAL_PATH="src/do3se_met"

# First check if repo is added
printf "====Attempting to remove old copies of the repositories"
git remote remove pyDO3SE
git remote remove thermal_time
git remote remove do3se_phenology
git remote remove do3se_met
printf "Removed old copies of the repositories=== \n\n"


# printf "====Attempting to add new copies of the repositories"
# git remote add -f pyDO3SE git@github.com:SEI-DO3SE/pyDO3SE
# git remote add -f thermal_time git@github.com:SEI-DO3SE/thermal_time
# git remote add -f do3se_phenology git@github.com:SEI-DO3SE/do3se_phenology
# git remote add -f do3se_met git@github.com:SEI-DO3SE/do3se_met
# printf "Added new copies of the repositories=== "


PYDO3SE_REMOTE="git@github.com:SEI-DO3SE/pyDO3SE"


git subtree add --prefix $PYDO3SE_LOCAL_PATH $PYDO3SE_REMOTE RELEASE --squash