import os
os.environ["NUMBA_DISABLE_JIT"] = "1"

from scripts.validators.common import VALIDATION_OUTPUT_DIR
from scripts.validators.figure_reproduction import (reproduce_2024_fig_16,
                                             reproduce_2025_fig_8)
from scripts.validators.octree_validation import generate_octree_figures


def validation_main() -> None:
    reproduce_2025_fig_8()
    reproduce_2024_fig_16()
    generate_octree_figures(n=5)


if __name__ == '__main__':
    VALIDATION_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    validation_main()
