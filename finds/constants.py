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
VALIDATION_GRAPH_PATH='v-2026-8.png'
VALIDATION_3D_GRAPH_PATH='v-2025-16.png'

DATA_FILE_NAME = 'data.h5'
ANIMATION_FILE_NAME = 'animation.mp4'
RADIAL_DENSITY_DISTRIBUTION_FILE_NAME = 'density_animation.mp4'
