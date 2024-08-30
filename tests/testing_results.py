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

scenarios = [
    "vanKeken1997_case1a",
    "Crameri2012_case2",
    "continental_rift",
    "punch",
    "Schmalholz_2011",
    "two_layers_with_variable_conductivity",
]

# Name of the files to compare
fields = [
    "time",
    "density",
    "heat",
    "pressure",
    "strain",
    "strain_rate",
    "surface",
    "temperature",
    "velocity",
    "viscosity",
    "step",
]

# steps to compare
steps = [0, 1, 10000]

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


@pytest.mark.parametrize("step", steps)
@pytest.mark.parametrize("field", fields)
@pytest.mark.parametrize("scenario", scenarios)
def test_result(scenario, field, step):
    """Run tests"""

    if scenario != 'Crameri2012_case2' and field == 'surface':
        pytest.skip('No surface for this scenario')

    if scenario == 'Crameri2012_case2' and step == 1:
        pytest.skip('Tested with only one processor')

    if (scenario == 'punch' and step == 1) or (scenario == 'punch' and field == 'step_1'):
        pytest.skip('Test with only one step')

    if (scenario != 'two_layers_with_variable_conductivity' and step == 10000) or (scenario == 'two_layers_with_variable_conductivity' and step == 1):
        pytest.skip('Test step/file not available for this scenario')

    test_path = base_path / "data" / scenario/ "output"
    expected_path = base_path / "data" / scenario/ "expected"

    filenames = [f"{field}_{step}.txt"]

    if field == 'step':
        if scenario in ['Crameri2012_case2']:
            filenames = [f"step_{step}_0.txt"]
        else:
            filenames = [f"step_{step}_0.txt", f"step_{step}_1.txt"]

    for filename in filenames:
        output = read(test_path / filename)
        expected = read(expected_path / filename)

        npt.assert_allclose(output, expected, rtol=2e-4, atol=1.0e-18)
