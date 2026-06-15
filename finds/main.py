from .fish import generate_system
from .postprocessing import process_data
from .simulation import perform_simulation
from tqdm import tqdm

import numpy as np

def main():
    system = generate_system(
        distribution='lattice',
        orientation='aligned',
        angle_delta=np.pi/2,
        size=6,
        spacing=5,
        debug_print=True
    )

    output_dir = perform_simulation(
        system,
        time_step=0.01,
        end_time=20,
        use_barnes_hut=True,
        bh_ratio=0.9,
        print_iterations=True,
        print_each_fish=False
    )

    process_data(output_dir)


if __name__ == '__main__':
    main()
