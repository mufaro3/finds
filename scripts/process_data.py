from pathlib import Path

import sys
import h5py
from matplotlib import rc
from tqdm import trange
import numpy as np

from processingmodules.TrajectoryPlot import TrajectoryPlot
from processingmodules.AnimationGenerator import AnimationGenerator
from processingmodules.CrossSectionGenerator import CrossSectionGenerator
from processingmodules.DensityAnimationGenerator import \
    DensityAnimationGenerator
from processingmodules.MeanRadialDistancePlot import MeanRadialDistancePlot
from processingmodules.ProcessingModule import ProcessingModule

DEFAULT_MODULES_LIST=[
    AnimationGenerator(),
    DensityAnimationGenerator(),
    CrossSectionGenerator(),
    MeanRadialDistancePlot(),
    TrajectoryPlot()
]

SYSTEM_DTYPE = np.dtype([
    ("positions", np.float64, (3,)),
    ("orientations", np.float64, (3,)),
    ("lengths", np.float64),
    ("sigmas", np.float64),
])

def process_data(
        output_dir: Path,
        data_file: Path,
        modules: list[ProcessingModule] = DEFAULT_MODULES_LIST,
        use_pdf: bool = False) -> None:
    r"""
    Processes the data stored within the output directory. Each state and
    simulation time are read in sequentially, one at a time, so the data is
    processed using "processing modules" that append one-at-a-time (to reduce
    the load on system memory).

    :param output_dir: The output directory to read from/write to.
    :type output_dir: Path

    :param generate_animation: Whether or not to generate an animation of the
      simulation.
    :type  generate_animation: bool

    :param use_pdf: Whether or not to output all image files to PDF (intended
      for inclusion in PDF documents, like reports). False by default.
    :type  use_pdf: bool
    """
    fontconfig = {
        'family': 'serif',
        'serif': ['Computer Modern'],
        'size': 14
    }
    rc('font', **fontconfig)
    rc('text', usetex=True)
    rc('svg', fonttype='none')

    with h5py.File(data_file, "r") as h5:
        time_dataset = h5["/time"]
        position_dataset = h5["/position"]
        orientation_dataset = h5["/orientation"]
        length_dataset = h5["/length"]
        sigma_dataset = h5["/volumetric_flow_rate"]

        num_frames = len(time_dataset)

        for module in modules:
            module.begin(output_dir, use_pdf)

        for frame_index in trange(num_frames, desc="Postprocessing"):
            time = time_dataset[frame_index]

            N = len(length_dataset[frame_index])
            system = np.empty(N, dtype=SYSTEM_DTYPE).view(np.recarray)

            system["positions"] = position_dataset[frame_index]
            system["orientations"] = orientation_dataset[frame_index]
            system["lengths"] = length_dataset[frame_index]
            system["sigmas"] = sigma_dataset[frame_index]

            for module in modules:
                module.append_state(
                    system,
                    time,
                    frame_index + 1,
                    num_frames
                )

        for module in modules:
            module.end()


if __name__ == '__main__':
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print(f"Usage: {sys.argv[0]} <data-file>")
        sys.exit(1)

    data_file = Path(sys.argv[1])

    if not data_file.exists():
        print(f"File does not exist: {data_file}")
        sys.exit(1)

    # Process into the same directory as the input file
    output_dir = data_file.parent

    process_data(output_dir, data_file)
