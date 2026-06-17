import shutil
from pathlib import Path
import argparse

from tqdm import tqdm, trange
import numpy as np
from matplotlib import pyplot as plt
from numpy.typing import NDArray
import time
from itertools import product, combinations

from finds.fish import generate_system
from finds.constants import DATA_FILE_NAME, VALIDATION_OUTPUT_PATH
from finds.io import close_filestream, init_input_filestream
from finds.simulation import perform_simulation
from finds.util import rejoin, split
from finds.calculations import calculate_system_derivative, OctreeNode, \
    build_octree

def generate_comparison_figure(output_dir: Path) -> None:
    """
    Generate timing comparison between brute-force and Barnes-Hut.
    """

    output_dir.mkdir(exist_ok=True, parents=True)

    def perform_time_test(n: int,
                          bh_ratio: float | None = None,
                          repeats: int = 3) -> float:
        """
        Time a single derivative calculation.
        """

        example_system = generate_system(
            distribution='random',
            orientation='random',
            n_random=n,
            debug_print=False
        )

        use_barnes_hut = bh_ratio is not None

        # warmup to exclude numba compilation
        calculate_system_derivative(
            example_system,
            use_barnes_hut,
            bh_ratio if bh_ratio is not None else 0.0
        )

        times = []

        for _ in range(repeats):
            start_time = time.perf_counter()

            calculate_system_derivative(
                example_system,
                use_barnes_hut,
                bh_ratio if bh_ratio is not None else 0.0
            )

            end_time = time.perf_counter()
            times.append(end_time - start_time)

        return np.mean(times)

    # system sizes
    lognvalues = np.flip(np.arange(1,3))

    # Barnes-Hut opening angles
    bh_ratios = np.arange(0.25, 1.0, step=0.25)

    fig, ax = plt.subplots(figsize=(10, 6))

    # brute-force baseline
    brute_force_times = []

    for logn in tqdm(lognvalues, desc='Brute Force Benchmark'):
        brute_force_times.append(
            perform_time_test(2 ** logn, bh_ratio=None)
        )

    ax.plot(
        lognvalues,
        brute_force_times,
        label="Brute force",
        linewidth=3
    )

    # Barnes-Hut curves
    for bh_ratio in tqdm(bh_ratios, leave=False, desc='Barnes-Hut'):
        bh_times = []

        for logn in tqdm(lognvalues, leave=False, desc=f'BH Ratio={bh_ratio:.2f}'):
            bh_times.append(perform_time_test(2 ** logn, bh_ratio=bh_ratio))

        ax.plot(
            lognvalues,
            bh_times,
            label=f"BH ratio={bh_ratio:.1f}"
        )

    ax.set_xlabel(r"Number of fish, $\log_2 N$")
    ax.set_ylabel("Runtime (seconds)")
    ax.set_title("Barnes–Hut vs Brute Force Runtime")
    ax.legend()
    ax.grid(True)

    output_path = output_dir / "barnes_hut_comparison.png"

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    print(f"Saved comparison figure to {output_path}")

if __name__ == '__main__':
    output_dir = Path(f'output/benchmark')
    output_dir.mkdir(parents=True, exist_ok=True)
    generate_comparison_figure(output_dir)
