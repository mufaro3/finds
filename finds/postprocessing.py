from pathlib import Path

from matplotlib import rc
from tqdm import trange

from finds.processingmodules.TrajectoryPlot import TrajectoryPlot

from .constants import DATA_FILE_NAME
from .io import close_filestream, init_input_filestream
from .processingmodules.AnimationGenerator import AnimationGenerator
from .processingmodules.CrossSectionGenerator import CrossSectionGenerator
from .processingmodules.DensityAnimationGenerator import \
    DensityAnimationGenerator
from .processingmodules.MeanRadialDistancePlot import MeanRadialDistancePlot
from .processingmodules.ProcessingModule import ProcessingModule

DEFAULT_MODULES_LIST=[
    AnimationGenerator(),
    DensityAnimationGenerator(),
    CrossSectionGenerator(),
    MeanRadialDistancePlot(),
    TrajectoryPlot()
]


def process_data(
        output_dir: Path,
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

    datafile = output_dir / DATA_FILE_NAME
    io = init_input_filestream(datafile, cache_limit=1024**2)

    for module in modules:
        module.begin(output_dir, use_pdf)

    num_frames = len(io.time_dataset)

    try:
        for frame_index in trange(num_frames, desc='Postprocessing'):
            system = io.state_dataset[frame_index]
            time = io.time_dataset[frame_index]

            for module in modules:
                module.append_state(system, time, frame_index + 1, num_frames)

    finally:
        for module in modules:
            module.end()

    close_filestream(io)
