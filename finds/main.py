from .fish import generate_system
from .postprocessing import process_data
from .simulation import perform_simulation


def main():
    system = generate_system(
        distribution='lattice',
        orientation='aligned'
    )

    output_dir = perform_simulation(
        system,
        time_step=0.1,
        end_time=50
    )

    process_data(
        output_dir,
        generate_animation=True,
        generate_density_animation=True
    )


if __name__ == '__main__':
    main()
