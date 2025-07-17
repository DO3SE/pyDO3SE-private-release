#!/bin/bash
set -e

source .venv/bin/activate

uv sync --all-groups

# Copy integration tests from pyDO3SE
rm -R tests/key_processes
mkdir -p tests/key_processes
touch tests/__init__.py
cp -R src/pyDO3SE/tests/key_processes tests
cp src/pyDO3SE/tests/utils.py tests/utils.py

rm -R examples
mkdir -p examples
cp -R src/pyDO3SE/examples .

# Delete all .pyc files
find tests -name "*.pyc" -type f -delete