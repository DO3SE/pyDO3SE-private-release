#!/bin/bash
set -e

uv venv --python=3.12
source .venv/bin/activate

uv sync --extra dev
uv sync --extra grid


# Copy integration tests from pyDO3SE
mkdir -p tests/key_processes
cp -R src/pyDO3SE/tests/key_processes tests

mkdir -p examples
cp -R src/pyDO3SE/examples .
