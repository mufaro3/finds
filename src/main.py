import numpy as np
from pprint import pprint
from .simulation import *
from .fish import *

FISH_LENGTH=1 # m
VOLUMETRIC_FLOW_RATE=4*np.pi # rad/s
FISH_ORIENTATION_DELTA=0.01 # rad
NUMBER_OF_FISH=3
FISH_GENERATION_BOUNDS=[10,10,10]
FISH_SELF_PROPELLED_SPEED=VOLUMETRIC_FLOW_RATE/(4*np.pi*FISH_LENGTH**2)

TIME_STEP=0.01
END_TIME=1
        
def main():
    system = generate_system(
        NUMBER_OF_FISH,
        FISH_GENERATION_BOUNDS,
        FISH_ORIENTATION_DELTA
    )
    pprint(system)
        
if __name__ == '__main__':
    main()
