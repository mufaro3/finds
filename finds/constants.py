from pathlib import Path

import numpy as np

#: Characteristic fish length (m)
FISH_LENGTH: float = 1.0

#: Volumetric flow rate (rad/s)
VOLUMETRIC_FLOW_RATE: float = 4 * np.pi

#: Self-propelled swimming speed (m/s)
FISH_SELF_PROPELLED_SPEED: float = (
    VOLUMETRIC_FLOW_RATE /
    (4 * np.pi * FISH_LENGTH**2)
)

#: The root output path
ROOT_OUTPUT_PATH=Path('output')

#: The output path for all of the simulation files
SIMULATION_OUTPUT_DIR=ROOT_OUTPUT_PATH / Path('simulations')

#: The output path for simulation files
SIMULATION_OUTPUT_NAME='sim'

#: The output path for validation files
VALIDATION_OUTPUT_PATH=ROOT_OUTPUT_PATH / Path('validation')

#: The output path for benchmark files
BENCHMARK_OUTPUT_PATH=ROOT_OUTPUT_PATH / Path('benchmark')

#: The name of all calculation data files
DATA_FILE_NAME='data.h5'
