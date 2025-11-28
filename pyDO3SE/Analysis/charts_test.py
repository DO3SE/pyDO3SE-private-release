import pytest
import pandas as pd
import numpy as np
from pyDO3SE.Output.Output_Shape import output_fields_map
from pyDO3SE.Analysis.charts import annual_graph, monthly_diurnal_graph


def test_create_an_annual_graph():
    # We can only test that this runs successfully.
    # If there are any actual differences these should appear in git diff
    fig = annual_graph([1, 2, 3, 8, 9, 7], output_fields_map["gsto"])
    assert fig


@pytest.mark.parametrize("month", [0, 6, 11])
def test_create_an_diurnal_graph(month):
    # We can only test that this runs successfully.
    # If there are any actual differences these should appear in git diff

    data = pd.DataFrame(
        {
            "hr": np.arange(0, 24),
            "gsto": np.sin(np.linspace(0, 2 * np.pi, 24)) * 5 + 5,
        }
    )
    observed_data = pd.DataFrame(
        {
            "hr": np.arange(0, 24),
            "observed": np.sin(np.linspace(0, 2 * np.pi, 24)) * 5 + 5 + np.random.normal(0, 1, 24),
        }
    )
    fig = monthly_diurnal_graph(
        data,
        observed_data,
        output_fields_map["gsto"],
        month,
        ylim=(0, 10),
    )
    assert fig is not None


if __name__ == "__main__":
    test_create_an_annual_graph()
