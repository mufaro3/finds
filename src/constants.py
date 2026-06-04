import numpy as np
from types import SimpleNamespace

Constants = SimpleNamespace(
    FISH_LENGTH=1, # m
    VOLUMETRIC_FLOW_RATE=4*np.pi, # rad/s
)

Constants.FISH_SELF_PROPELLED_SPEED=\
    Constants.VOLUMETRIC_FLOW_RATE/(4*np.pi*Constants.FISH_LENGTH**2)
