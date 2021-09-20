"""
Compare the Mandyoc output data with the expected data using the Crameri model 
"""
import os
import pytest
import numpy as np
import numpy.testing as npt
from pathlib import Path

# Test path
base_path = Path(os.path.realpath(os.path.abspath(__file__))).parent
test_path = base_path / "testing_data"

# Expected data
# expected_path = os.path.join(test_path, "expected")
expected_path = test_path / "expected"

# Name of the files to compare
fields = [
    "density",
    "heat",
    "pressure",
    "strain",
    "strain_rate",
    "temperature",
    "sp_surface_global",
    "viscosity",
    "step_0",
    "step_1",
    "time",
    "velocity",
]
# steps to compare
steps = [0, 1]


def read(filename):
    """
    Read the file
    """
    if "time" not in filename.name:
        args = dict(unpack=True, comments="P", skiprows=2)
    else:
        args = dict(unpack=True, delimiter=":", usecols=1)
    data = np.loadtxt(filename, **args)
    return data


@pytest.mark.parametrize("field", fields)
@pytest.mark.parametrize("step", steps)
def test_result(field, step):
    """ """
    filename = f"{field}_{step}" + ".txt"
    output = read(test_path / filename)
    expected = read(expected_path / filename)

    npt.assert_allclose(output, expected, rtol=1e-5)
