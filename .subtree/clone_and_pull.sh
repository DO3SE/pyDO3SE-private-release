# Initial clone
PYDO3SE_LOCAL_PATH="src/pyDO3SE"
THERMAL_TIME_LOCAL_PATH="src/thermal_time"
DO3SE_PHENOLOGY_LOCAL_PATH="src/do3se_phenology"
DO3SE_MET_LOCAL_PATH="src/do3se_met"

# First check if repo is added
git remote remove pyDO3SE
git remote remove thermal_time
git remote remove do3se_phenology
git remote remove do3se_met

git remote add -f pyDO3SE git@github.com/SEI-DO3SE/pyDO3SE
git remote add -f thermal_time git@github.com/SEI-DO3SE/thermal_time
git remote add -f do3se_phenology git@github.com/SEI-DO3SE/do3se_phenology
git remote add -f do3se_met git@github.com/SEI-DO3SE/do3se_met

# TODO: Might need to use release branch
# TODO: Fix if statement
if [$1 = "clone"]; then
  git subtree add --prefix $PYDO3SE_LOCAL_PATH pyDO3SE main --squash
  git subtree add --prefix $THERMAL_TIME_LOCAL_PATH thermal_time main --squash
  git subtree add --prefix $DO3SE_PHENOLOGY_LOCAL_PATH do3se_phenology main --squash
  git subtree add --prefix $DO3SE_MET_LOCAL_PATH do3se_met main --squash
fi
if [$1 == "pull"]; then
  # Update subtree
  git subtree pull --prefix $PYDO3SE_LOCAL_PATH pyDO3SE main --squash
  git subtree pull --prefix $THERMAL_TIME_LOCAL_PATH thermal_time main --squash
  git subtree pull --prefix $DO3SE_PHENOLOGY_LOCAL_PATH do3se_phenology main --squash
  git subtree pull --prefix $DO3SE_MET_LOCAL_PATH do3se_met main --squash
fi
if [$1 == "tag"]; then
# TODO: Update version info
# Should store the open version somewhere

# We can either create a specific hash of the versions of each dependency or just use the pyDO3SE version


# == User Documentation == #
# TODO: Where should we store the end user documentation? This probably needs to combine documentation from all repositories
