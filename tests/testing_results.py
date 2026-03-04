"""
Compare the Mandyoc output data with the expected data
"""
import pytest
import numpy as np
import numpy.testing as npt
from pathlib import Path


BASE_PATH = Path(__file__).parent.absolute()

MAIN_FIELDS = [
    "density",
    "heat",
    "pressure",
    "step",
    "strain",
    "strain_rate",
    "temperature",
    "time",
    "velocity",
    "viscosity",
]


# Definition of each scenario
SCENARIOS = {
    "vanKeken1997_case1a": {
        "steps": [0, 1],
        "fields": MAIN_FIELDS,
    },
    "Crameri2012_case2": {
        "steps": [0],
        "fields": MAIN_FIELDS,
        "num_procs": 1,
    },
    "continental_rift": {
        "steps": [0, 1],
        "fields": MAIN_FIELDS,
    },
    "punch": {
        "steps": [0],
        "fields": MAIN_FIELDS,
    },
    "Schmalholz_2011": {
        "steps": [0, 1],
        "fields": MAIN_FIELDS,
    },
    "two_layers_with_variable_conductivity": {
        "steps": [0, 10000],
        "fields": MAIN_FIELDS,
    },
}


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


def get_test_cases():
    """Generates valid test combinations."""
    cases = []
    for name, config in SCENARIOS.items():
        for step in config["steps"]:
            for field in config["fields"]:
                cases.append((name, field, step))
    return cases


@pytest.mark.parametrize("scenario, field, step", get_test_cases())
def test_result(scenario, field, step):
    config = SCENARIOS[scenario]

    test_path = BASE_PATH / "data" / scenario / "output"
    expected_path = BASE_PATH / "data" / scenario / "expected"

    if not test_path.exists():
        pytest.fail(f"Simulation directory missing. Scenario {scenario} likely failed to run.")

    # Handle the 'step' field
    if field == 'step':
        num_procs = config.get("num_procs", 2)
        filenames = [f"step_{step}_{i}.txt" for i in range(num_procs)]
    else:
        filenames = [f"{field}_{step}.txt"]

    for fname in filenames:
        out_file = test_path / fname
        exp_file = expected_path / fname

        assert out_file.exists(), f"Missing output file: {out_file}"

        output = read(out_file)
        expected = read(exp_file)

        npt.assert_allclose(output, expected, rtol=2e-4, atol=1.0e-18)
