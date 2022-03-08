"""
Compare the Mandyoc output data with the expected data
"""
import os
import pytest
import numpy as np
import numpy.testing as npt
from pathlib import Path

# Test path
base_path = Path(os.path.realpath(os.path.abspath(__file__))).parent

# Name of the files to compare
fields = [
    "time",
    "density",
    "heat",
    "pressure",
    "strain",
    "strain_rate",
    "temperature",
    "velocity",
    "viscosity",
    "step_0",
    "step_1",
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
    """Run tests"""

    test_path = base_path / "output"
    expected_path = base_path / "expected"

    filename = f"{field}_{step}" + ".txt"
    output = read(test_path / filename)
    expected = read(expected_path / filename)

    npt.assert_allclose(output, expected, rtol=2e-4)
