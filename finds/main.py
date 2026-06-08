from pprint import pprint

from .fish import generate_system
from .simulation import perform_simulation 

def main():
    system = generate_system(
        n=3,
        bounds=[10,10,10],
        angle_delta=0.01
    )

    output_file = perform_simulation(
        system,
        time_step=0.1,
        end_time=1
    )

if __name__ == '__main__':
    main()
