from .fish import generate_system
from .postprocessing import process_data
from .simulation import perform_simulation


def main():
    print('Generating system..')
    system = generate_system(
        distribution='lattice',
        orientation='aligned',
        size=20,
        spacing=5.0
    )
    print(f'System size: N={system.shape[0]}')

    output_dir = perform_simulation(
        system,
        time_step=0.1,
        end_time=50
    )

    print('Generating figures..')
    process_data(
        output_dir,
        generate_animation=True,
        generate_density_animation=True
    )


if __name__ == '__main__':
    main()
