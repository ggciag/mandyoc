"""
Testing Mandyoc code using Crameri example for only 3 time steps
"""

import numpy as np
import os

# Test path
test_path = "test/testing_data"
# Path to the expected data to make the comparison
expected_result_path = os.path.join(test_path, "expected")
# Name of the files to compare
file_name = [
    "density_{}",
    "heat_{}",
    "pressure_{}",
    "strain_{}",
    "strain_rate_{}",
    "temperature_{}",
    #"sp_surface_global_{}",
    "viscosity_{}",
    "step_{}_0",
    "step_{}_1",
    "time_{}",
    "velocity_{}",
]
# steps to compare
steps = [0, 1]

rtol = 0.1 

# Error message
error_message = "The test result is not equal or similar to the expected value for {}. Please contact the Mandyoc developers."

for step in steps:
    for name in file_name:
        if name != "time_{}":
            # Load the Mandyoc result
            filename = os.path.join(test_path, name + ".txt").format(step)
            data_step = np.loadtxt(
                filename,
                unpack=True,
                comments="P",
                skiprows=2,
            )

            # Load the expected result
            expected_filename = os.path.join(
                expected_result_path, name + ".txt"
            ).format(step)
            expected_data_step = np.loadtxt(
                expected_filename,
                unpack=True,
                comments="P",
                skiprows=2,
            )
        else:
            # Load the Mandyoc result for time
            filename = os.path.join(test_path, name + ".txt").format(step)
            data_step = time = np.loadtxt(
                filename, unpack=True, delimiter=":", usecols=(1)
            )
            # data_step +=1
            # Load the expected result for time
            expected_filename = os.path.join(
                expected_result_path, name + ".txt"
            ).format(step)
            expected_data_step = time = np.loadtxt(
                expected_filename, unpack=True, delimiter=":", usecols=(1)
            )

        # Compare the Mandyoc results with the expected results
        if (abs(data_step - expected_data_step) >= rtol).any():
            raise ValueError(error_message.format(filename))

print("Mandyoc test has run successfully.")
