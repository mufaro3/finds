from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import override
from tqdm import trange, tqdm

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FFMpegWriter
from numpy.typing import ArrayLike, NDArray

from .constants import DATA_FILE_NAME
from .io import close_filestream, init_input_filestream
from .util import split


class ProcessingModule(ABC):
    r"""
    This is an abstract class for all "processing modules" so that they can
    operate in parallel to reduce the amount of data read to the disk at one
    time (such as generating figures, performing calculations, etc).
    """
    output_dir: Path

    def begin(self, output_dir: Path):
        """
        Instructs the processing module to initialize such that it can
        produce data.

        :param output_dir: The output directory to write to.
        :type  output_dir: Path
        """
        self.output_dir = output_dir

    @abstractmethod
    def append_state(self, system: NDArray, time: float) -> None:
        r"""
        Appends the state :math:`\mathbf{X}` and time :math:`t` to the current
        process being generated (whether by animating it, including it in a
        computation, etc.).

        :param system: The state of the system at time :math:`t`.
        :type  system: NDArray

        :param time: The time, :math:`t`.
        :type  time: float
        """
        pass

    def end(self) -> None:
        pass


@dataclass
class AnimationGenerator(ProcessingModule):
    r"""
    Generates a 3-D animation of the full fish simulation, taking (0,0) to be
    the origin. Allows for rapid customization of simulation colors, viewing,
    and debug parameters as well.
    """

    # sizing and configuration
    particle_radius: int = 30
    orientation_width: int = 5
    orientation_length: int = 5
    padding: int = 0.1
    figsize: tuple[int] = (8,8)

    # color scheme
    head_color: str = 'tab:red'
    tail_color: str = 'tab:blue'
    particle_color: str = 'tab:cyan'
    orientation_color: str = 'tab:orange'

    def __init__(self,
                 fps: int = 60,
                 max_bounds: ArrayLike = [ 100, 100, 100 ]):
        self.fps = fps
        self.max_bounds = np.asarray(max_bounds, dtype=float)

    def set_limits(self):
        xmin, ymin, zmin = -self.max_bounds
        xmax, ymax, zmax =  self.max_bounds

        self.ax.set_xlim(xmin - self.padding, xmax + self.padding)
        self.ax.set_ylim(ymin - self.padding, ymax + self.padding)
        self.ax.set_zlim(zmin - self.padding, zmax + self.padding)

        self.ax.set_box_aspect([
            xmax - xmin,
            ymax - ymin,
            zmax - zmin
        ])

    @override
    def begin(self, output_dir: Path) -> None:
        r"""
        Sets up the MatPlotLib animation for rendering.

        :type max_bounds: ArrayLike
        """
        super().begin(output_dir)

        self.fig = plt.figure(figsize=self.figsize)
        self.ax = self.fig.add_subplot(projection='3d')
        self.set_limits()
        self.ax.view_init(elev=20, azim=45)

        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        self.ax.set_zlabel("Z")

        # initialize the position particles
        self.scatter_plot = self.ax.scatter(
            [], [], [],
            s=self.particle_radius,
            c=self.particle_color)

        # initialize the orientations
        self.quiver_plot = self.ax.quiver(
            0, 0, 0, 0, 0, 0,
            color=self.orientation_color)

        self.title_text = self.ax.set_title('')
        self.output_filename = self.output_dir / '3d_animation.mp4'
        self.writer = FFMpegWriter(self.fps)
        self.writer.setup(self.fig, str(self.output_filename), dpi=100)

    @override
    def append_state(self, system: NDArray, time: float) -> None:
        r"""
        Draws the current state to the animation.
        """
        positions, orientations = split(system)

        self.scatter_plot._offsets3d = (
            positions[:, 0],
            positions[:, 1],
            positions[:, 2]
        )

        vecs = orientations * self.orientation_length
        self.quiver_plot.set_segments(
            np.stack([positions, positions + vecs], axis=1)
        )

        self.title_text.set_text(f't={time:.2f}')
        self.writer.grab_frame()

    @override
    def end(self) -> None:
        r"""
        Closes the animation and saves the file.
        """
        self.writer.finish()
        plt.close(self.fig)
        print(f'Animation saved at {self.output_filename}')


@dataclass
class DensityAnimationGenerator(ProcessingModule):
    r"""
    Produces an animation of the radial mass distribution of fish from the
    center-of-mass for each time step, :math:`\delta t`.
    """
    n_bins: int = 30
    max_radius: float = 50.0  # adjust or auto-compute if you prefer

    @override
    def __init__(self, fps: int = 60):
        self.fps = fps

    @override
    def begin(self, output_dir: Path) -> None:
        self.output_dir = output_dir
        self.fig, self.ax = plt.subplots(figsize=(6, 4))
        self.bins = np.linspace(0, self.max_radius, self.n_bins + 1)
        bin_centers = 0.5 * (self.bins[:-1] + self.bins[1:])

        self.bars = self.ax.bar(
            bin_centers,
            np.zeros(self.n_bins),
            width=self.bins[1] - self.bins[0])

        self.ax.set_xlim(0, self.max_radius)
        self.ax.set_xlabel("Radius")
        self.ax.set_ylabel("Count")
        self.title_text = self.ax.set_title("Radial Density Distribution")

        self.output_filename = self.output_dir /\
            'radial_density_evolution.mp4'
        self.writer = plt.matplotlib.animation.FFMpegWriter(fps=self.fps)
        self.writer.setup(self.fig, str(self.output_filename), dpi=100)

    @override
    def append_state(self, system: NDArray, time: float) -> None:
        positions, _ = split(system)
        com = np.mean(positions, axis=0)
        radii = np.linalg.norm(positions - com, axis=1)
        counts, _ = np.histogram(radii, bins=self.bins)

        for bar, count in zip(self.bars, counts):
            bar.set_height(count)

        current_ymax = self.ax.get_ylim()
        max_count = np.max(counts)
        if np.any(max_count > current_ymax):
            self.ax.set_ylim(0, max_count * 1.1)

        self.title_text.set_text(f"Radial Density (t={time:.2f})")
        self.writer.grab_frame()

    @override
    def end(self) -> None:
        """
        Closes the animation and saves the file.
        """
        if self.writer is not None:
            self.writer.finish()
        plt.close(self.fig)

        print('Saved density distribution animation to '+\
              f'{self.output_filename}')


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
        total_iterations = len(io.time_dataset)
        interval = max(total_iterations // 10, 1)

        for i in trange(total_iterations, desc='Postprocessing'):
            system = io.state_dataset[i]
            time = io.time_dataset[i]

            for module in modules:
                module.append_state(system, time)

    finally:
        for module in modules:
            module.end()

    close_filestream(io)
