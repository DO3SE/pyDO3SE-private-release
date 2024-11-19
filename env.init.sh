uv venv --python=3.12
source .venv/bin/activate
uv pip install -r requirements/common.txt
uv pip install -r requirements/tests.txt
uv pip install -r requirements/local.txt --no-deps # These should be installed separately to ensure they do not clash with the other requirements
