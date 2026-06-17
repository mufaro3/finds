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

#: The output path for validation files
VALIDATION_OUTPUT_PATH='validation'

#: The name of all calculation data files
DATA_FILE_NAME = 'data.h5'
