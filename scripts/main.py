from finds.fish import generate_system
from finds.postprocessing import process_data
from finds.simulation import perform_simulation


def main():
    system = generate_system(
        distribution='lattice',
        orientation='swirl',
        angle_delta=0,
        size=100,
        spacing=25,
        debug_print=True
    )

    output_dir = perform_simulation(
        system,
        end_time=10,

        # Barnes-Hut
        use_barnes_hut=True,
        bh_ratio=0.5,

        # Integration
        integration_method = 'RK45',
        rtol = 1e-6,
        atol = 1e-8,
        time_step = 0.01,

        # Debug Printing
        print_time_progression=True,
        print_each_fish=False
    )

    process_data(output_dir)


if __name__ == '__main__':
    main()
