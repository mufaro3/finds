from pathlib import Path

from tqdm import trange

from .constants import DATA_FILE_NAME
from .io import close_filestream, init_input_filestream
from .processingmodules.AnimationGenerator import AnimationGenerator
from .processingmodules.DensityAnimationGenerator import \
    DensityAnimationGenerator
from .processingmodules.ProcessingModule import ProcessingModule

DEFAULT_MODULES_LIST=[
    AnimationGenerator(),
    DensityAnimationGenerator()
]


def process_data(
        output_dir: Path,
        modules: list[ProcessingModule] = DEFAULT_MODULES_LIST) -> None:
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
    """
    datafile = output_dir / DATA_FILE_NAME
    io = init_input_filestream(datafile, cache_limit=1024**2)

    for module in modules:
        module.begin(output_dir)

    try:
        for i in trange(len(io.time_dataset), desc='Postprocessing'):
            system = io.state_dataset[i]
            time = io.time_dataset[i]

            for module in modules:
                module.append_state(system, time)

    finally:
        for module in modules:
            module.end()

    close_filestream(io)
