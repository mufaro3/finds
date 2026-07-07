from pathlib import Path

import numpy as np

#: The root output path
ROOT_OUTPUT_PATH=Path('output')

#: The output path for all of the simulation files
SIMULATION_OUTPUT_DIR=ROOT_OUTPUT_PATH / Path('simulations')

#: The output directory for the latest simulation
SIMULATION_LATEST_DIR='latest_simulation'

#: The output path for simulation files
SIMULATION_OUTPUT_NAME='sim'

#: The output path for validation files
VALIDATION_OUTPUT_PATH=ROOT_OUTPUT_PATH / Path('validation')

#: The output path for benchmark files
BENCHMARK_OUTPUT_PATH=ROOT_OUTPUT_PATH / Path('benchmark')

#: The name of all calculation data files
DATA_FILE_NAME='data.h5'
