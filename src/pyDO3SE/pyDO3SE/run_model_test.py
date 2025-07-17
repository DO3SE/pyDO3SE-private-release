"""Unit tests for run_model.py ."""
import pytest

from pytest_mock import MockerFixture
from unittest.mock import MagicMock, Mock

from .run_model import run_model


__module_loc__ = __package__ + '.run_model'


def Obj(name: str, *args, **kwargs) -> object:
    argsAsKwargs = {k._extract_mock_name(): k for k in args}
    return type(name, (object,), {**kwargs, **argsAsKwargs})()


# == Mocked return values == #
mock_config = MagicMock()
mock_initial_state = MagicMock()
mock_external_state = MagicMock()
mock_model_processes = MagicMock()
mock_state_out = MagicMock()
mock_run_processes = MagicMock(return_value=mock_state_out)
mock_state_logs = [1, 2, 3]
# mock_processes = [1, 2, 3]
# mock_overrides = MagicMock()
mock_process_runner = Mock(
    run_processes=mock_run_processes,
    state_logs=mock_state_logs,
)


@pytest.fixture(autouse=True)
def mock_imports(mocker: MockerFixture) -> object:
    """Mock imported functions."""
    ProcessRunnerMock = mocker.patch(__module_loc__ + '.ProcessRunner',
                                     return_value=mock_process_runner)
    return Obj('MockedImports', ProcessRunnerMock)


class TestRunModel:
    @pytest.fixture(autouse=True)
    def _setup(self, mocker) -> None:
        pass

    def test_creates_a_process_runner_with_config_input(self, mock_imports):
        mock_config.output.fields = [1, 2, 3]
        out = run_model(mock_initial_state, mock_config, mock_external_state, mock_model_processes)
        mock_imports.ProcessRunner.assert_called_once_with(mock_config, mock_external_state, DEBUG_MODE=False)
        assert mock_process_runner.run_processes.call_count == 1
        assert out[0] == mock_state_out
        assert out[1] == mock_state_logs
