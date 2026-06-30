from finds.fish import generate_system
from finds.postprocessing import process_data
from finds.simulation import perform_simulation

def main():
    system = generate_system(
        distribution='sphere',
        orientation='swirl',
        angle_delta=0,
        size=50,
        spacing=10,
        debug_print=True
    )

    output_dir = perform_simulation(
        system,
        end_time=10,

        # Barnes-Hut
        use_barnes_hut=True,
        bh_ratio=0.75,

        # Integration
        integration_method = 'RK45',
        rtol = 1e-3,
        atol = 1e-4,
        time_step = 0.01,

        # Debug Printing
        print_time_progression=True,
        print_each_fish=False
    )

    process_data(output_dir)


if __name__ == '__main__':
    main()
