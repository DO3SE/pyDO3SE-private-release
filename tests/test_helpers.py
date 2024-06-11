import os
import pytest
from functools import wraps


def long_test(fn):
    @wraps(fn)
    def check_run_long(*args, **kwargs):
        if os.environ.get('TQUICK', 'False') == 'True':
            pytest.skip('Skip test in TQUICK mode as it takes a while to run')
        else:
            return fn(*args, **kwargs)
    return check_run_long
