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

VALIDATION_OUTPUT_PATH='validation'
DATA_FILE_NAME = 'data.h5'
