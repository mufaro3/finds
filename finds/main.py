from .fish import generate_system
from .simulation import perform_simulation
from .postprocessing import process_data


def main():
    system = generate_system(
        n=30,
        bounds=[1,1,0],
        angle_delta=0.01
    )

    output_dir = perform_simulation(
        system,
        time_step=0.1,
        end_time=100
    )

    process_data(
        output_dir,
        generate_animation=True,
        generate_density_animation=False
    )


if __name__ == '__main__':
    main()
