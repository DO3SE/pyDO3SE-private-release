
import warnings
from pyDO3SE.Config.config_validation import validate_config_fields
from data_helpers.list_helpers import flatten_list
from pyDO3SE.Config.Config_Shape import Config_Shape
from typing import List
from proflow.Objects.Process import Process
from proflow.Objects.Interface import I


def assert_greater_than(amount, strict=True):
    """Closure to assert inner less than amount."""
    def _assert_greater_than(**args):
        value = list(args.values())[0]
        name = list(args.keys())[0]
        try:
            assert value is not None, f"{name} is not defined"
            assert min(
                value) >= amount, f"Failed to assert {name} is greater than: {amount}. Min value is {min(value)}"
        except TypeError:
            raise TypeError(f"{name} is not a iterable")
        except Exception as e:
            if strict:
                raise e
            warnings.warn(str(e))
    return _assert_greater_than


def assert_less_than(amount):
    """Closure to assert inner less than amount."""
    def _assert_less_than(val):
        assert val < amount
    return _assert_less_than


def assert_list_defined(**args):
    value = list(args.values())[0]
    name = list(args.keys())[0]
    try:
        assert value is not None
        assert len(value) > 0
        assert value[0] is not None
    except AssertionError:
        raise ValueError(f"{name} is not defined")
    except TypeError:
        raise ValueError(f"{name} is not a iterable")
    except Exception as e:
        print(f"Failed to assert list: {name}")
        raise e


def validate_external_state(config: Config_Shape) -> List[Process]:
    assert_defined_processes = [
        Process(
            func=assert_list_defined,
            external_state_inputs=lambda e_state, row_index: [
                I(e_state.sinB, as_="sinB")]
        ),
        Process(
            func=assert_list_defined,
            external_state_inputs=lambda e_state, row_index: [
                I(e_state.VPD_dd, as_="VPD_dd")]
        ),
        Process(
            func=assert_list_defined,
            external_state_inputs=lambda e_state, row_index: [
                I(e_state.VPD, as_="VPD")]
        ),
        Process(
            func=assert_list_defined,
            external_state_inputs=lambda e_state, row_index: [
                I(e_state.PAR, as_="PAR")]
        ),
        Process(
            func=assert_list_defined,
            external_state_inputs=lambda e_state, row_index: [
                I(e_state.R, as_="R")]
        ),
        Process(
            func=assert_list_defined,
            external_state_inputs=lambda e_state, row_index: [
                I(e_state.O3, as_="O3")]
        ),
        Process(
            func=assert_greater_than(0, strict=False),
            external_state_inputs=lambda e_state, row_index: [
                I(e_state.VPD, as_="VPD")]
        ),
    ]

    return assert_defined_processes


def validate_config(config: Config_Shape) -> List[Process]:
    return [
        Process(
            func=validate_config_fields, config_inputs=lambda config: [
                I(config, as_='config'),
            ],
        ),
    ]


def validate_initial_state(config: Config_Shape) -> List[Process]:
    return []


def setup_validation_processes(config: Config_Shape) -> List[Process]:
    """Validates, config, state and exteranl data before running model"""
    return flatten_list([
        validate_config(config),
        validate_external_state(config),
        validate_initial_state(config),
    ])
