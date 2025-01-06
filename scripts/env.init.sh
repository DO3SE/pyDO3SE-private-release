#!/bin/bash
set -e

uv venv --python=3.12
source .venv/bin/activate
uv pip install -r requirements/common.txt
uv pip install -r requirements/tests.txt
uv pip install -r requirements/local.txt --no-deps # These should be installed separately to ensure they do not clash with the other requirements

# Copy integration tests from pyDO3SE
rm -R tests/key_processes
cp -R src/pyDO3SE/tests/key_processes tests/key_processes

rm -R examples
cp -R src/pyDO3SE/examples examples