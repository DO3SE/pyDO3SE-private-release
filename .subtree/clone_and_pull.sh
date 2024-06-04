# Initial clone
PYDO3SE_LOCAL_PATH="src/pyDO3SE"
THERMAL_TIME_LOCAL_PATH="src/thermal_time"
DO3SE_PHENOLOGY_LOCAL_PATH="src/do3se_phenology"
DO3SE_MET_LOCAL_PATH="src/do3se_met"

PYDO3SE_REMOTE="git@github.com:SEI-DO3SE/pyDO3SE"
THERMAL_TIME_REMOTE="git@github.com:SEI-DO3SE/thermal_time"
DO3SE_PHENOLOGY_REMOTE="git@github.com:SEI-DO3SE/do3se_phenology"
DO3SE_MET_REMOTE="git@github.com:SEI-DO3SE/do3se_met"

echo $1
echo $2


BRANCH=$2
# TODO: Might need to use release branch
# TODO: Fix if statement
if [[ "$1" == "clone" ]]; then
  echo "Cloning"
  git subtree add --prefix $PYDO3SE_LOCAL_PATH $PYDO3SE_REMOTE $BRANCH --squash
  git subtree add --prefix $THERMAL_TIME_LOCAL_PATH $THERMAL_TIME_REMOTE $BRANCH --squash
  git subtree add --prefix $DO3SE_PHENOLOGY_LOCAL_PATH $DO3SE_PHENOLOGY_REMOTE $BRANCH --squash
  git subtree add --prefix $DO3SE_MET_LOCAL_PATH $DO3SE_MET_REMOTE $BRANCH --squash
fi
if [[ "$1" == "pull" ]]; then
  echo "Pulling"
  # Update subtree
  git subtree pull --prefix $PYDO3SE_LOCAL_PATH $PYDO3SE_REMOTE $BRANCH --squash
  git subtree pull --prefix $THERMAL_TIME_LOCAL_PATH $THERMAL_TIME_REMOTE $BRANCH --squash
  git subtree pull --prefix $DO3SE_PHENOLOGY_LOCAL_PATH $DO3SE_PHENOLOGY_REMOTE $BRANCH --squash
  git subtree pull --prefix $DO3SE_MET_LOCAL_PATH $DO3SE_MET_REMOTE $BRANCH --squash
fi
if [[ "$1" == "tag" ]]; then
  # TODO: Complete tag implementation
  echo "Tagging"
  echo "Not implemented yet"
fi

# TODO: Update version info
# Should store the open version somewhere

# We can either create a specific hash of the versions of each dependency or just use the pyDO3SE version


# == User Documentation == #
# TODO: Where should we store the end user documentation? This probably needs to combine documentation from all repositories
