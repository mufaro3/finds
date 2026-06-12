from .fish import generate_system
from .postprocessing import process_data
from .simulation import perform_simulation


def main():
    print('Generating system..')
    system = generate_system(
        distribution='random',
        orientation='aligned',
        n_random=200,
        bounds=[10, 10, 10]
    )
    print(f'System size: N={system.shape[0]}')

    output_dir = perform_simulation(
        system,
        time_step=0.1,
        end_time=10,
        use_barnes_hut=True,
        bh_ratio=0.5
    )

    print('Generating figures..')
    process_data(
        output_dir,
        generate_animation=True,
        generate_density_animation=True
    )


if __name__ == '__main__':
    main()
