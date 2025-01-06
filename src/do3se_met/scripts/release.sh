[ -z "$1" ] && echo "make sure to input [patch|minor|major]" && exit 1;
source venv/bin/activate
if [ -z "$(git status --porcelain)" ]; then
  echo "Working directory clean"
else
    git stash save -u "hold_build_deploy"
fi
bumpversion $1

python -m build

git branch -d RELEASE
git checkout --orphan RELEASE
git commit -m "RELEASE"
git push --force --set-upstream origin RELEASE
git push --tags
git checkout main

if [ -z "$(git status --porcelain)" ]; then
  echo "Working directory clean"
else
    git stash pop
fi