from finds.fish import generate_system
from finds.postprocessing import process_data
from finds.simulation import perform_simulation


def main():
    system = generate_system(
        distribution='lattice',
        orientation='radial inward',
        angle_delta=0,
        size=100,
        spacing=25,
        debug_print=True
    )

    output_dir = perform_simulation(
        system,
        time_step=0.01,
        end_time=5,
        use_barnes_hut=False,
        bh_ratio=0.5,
        print_iterations=True,
        print_each_fish=False
    )

    process_data(output_dir)


if __name__ == '__main__':
    main()
