from typing import Any, Callable, Dict, List, Tuple
from unittest.mock import MagicMock
from unittest import mock as umock
from enum import Enum
import decorator
import json
from data_helpers.encoders import AdvancedJsonEncoder

from pytest_mock import MockerFixture

MockedFunctions = Dict[str, MagicMock]



class ETestTypes(Enum):

    UNIT_TEST = "Unit Test"
    INTEGRATION_TEST = "Integration Test"


def setup_test(test_type: ETestTypes):
    def _inner(func):
        def wrapper(func, self, *args, **kwargs):
            print(f"Test is [{test_type.value}]")
            if test_type == ETestTypes.UNIT_TEST:
                self.setup_for_unit_test()
                self.mock_for_unit_test(umock)
                out = func(self, *args, **kwargs)
                mock_funcs_teardown(self)
                return out
            if test_type == ETestTypes.INTEGRATION_TEST:
                self.setup_for_integration_test()
                self.mock_for_integration_test(umock)
                out = func(self, *args, **kwargs)
                mock_funcs_teardown(self)
                return out
            return func(*args, **kwargs)
        return decorator.decorator(wrapper, func)
    return _inner



def mock_funcs(
    self,
    mocker,
    module_loc: str,
    funcs: List[Tuple[str, Any, Callable[[Any], Any]]],
) -> List[MockedFunctions]:
    """Mock and spy on functions.

    This will add each function as mock_functionname to self.

    Parameters
    ----------
    mocker : pytest mocker
        pytest mocker
    module_loc : str
        module location of functions to mock
    funcs : List[Tuple[str, Any, Callable[[Any], Any]]]
        functions to mock as list of Tuples
        Each Tuple is the name of the function, the return value then optionally the side effect

    Returns
    -------
    List[MockedFunctions]
        returns mocked functions

    """
    mocked_funcs = {
        'mock_' +
        k: mocker.patch(
            module_loc + '.' + k,
            return_value=v
        ) if not s
        else mocker.patch(
            module_loc + '.' + k,
            side_effect=s[0]
        )
        for k, v, *s in funcs
    }

    for k, v in mocked_funcs.items():
        setattr(self, k, v.start())
    setattr(self, 'mocked_funcs', mocked_funcs)


def mock_funcs_teardown(self):
    for k, v in self.mocked_funcs.items():
        v.stop()

def process_snapshot(data):
    return json.dumps(data, cls=AdvancedJsonEncoder, indent=4, sort_keys=True)

