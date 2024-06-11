import pytest
from pytest_mock import MockerFixture

from pyDO3SE.Output.Output_Shape import output_fields_map
from pyDO3SE.Analysis.charts import annual_graph

__module_loc__ = __package__ + '.charts'


@pytest.fixture(autouse=True)
def mock_pyplot(mocker: MockerFixture) -> object:
    """Mock imported functions."""
    _mock_pyplot = mocker.patch(__module_loc__ + '.plt', autospec=True)
    return _mock_pyplot


def test_create_an_annual_graph_of_gsto(mock_pyplot):
    # We can only test that this runs successfully.
    # If there are any actual differences these should appear in git diff
    fig = annual_graph([1, 2, 3, 8, 9, 7], output_fields_map['gsto'])
    assert fig
    mock_pyplot.figure.assert_called_with()
    # assert fig.texts == [Text(0.5, 0.98, 'Gsto')]


if __name__ == "__main__":
    test_create_an_annual_graph_of_gsto()
