#!/bin/bash
set -e

MAIN_BRANCH="main"

[ -z "$1" ] && echo "make sure to input [patch|minor|major]" && exit 1;
source venv/bin/activate
if [ -z "$(git status --porcelain)" ]; then
  echo "Working directory clean"
else
    git stash save -u "hold_build_deploy"
fi
bumpversion $1

git branch -D RELEASE
git checkout --orphan RELEASE
git commit -m "RELEASE"
git push --force --set-upstream origin RELEASE
git checkout $MAIN_BRANCH

if [ -z "$(git status --porcelain)" ]; then
  echo "Working directory clean"
else
    git stash pop
fi