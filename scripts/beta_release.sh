source venv/bin/activate

if [ -z "$(git status --porcelain)" ]; then
  echo "Working directory clean"
else
    git stash save -u "hold_build_deploy"
fi

# TODO: Add beta tag to bump version
# bumpversion $1
RELEASE_BRANCH="BETA"
echo $RELEASE_BRANCH

git checkout --orphan $RELEASE_BRANCH
git push --force --set-upstream origin $RELEASE_BRANCH

if [ -z "$(git status --porcelain)" ]; then
  echo "Working directory clean"
else
    git stash pop
fi