from abc import ABC, abstractmethod
from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.animation import FFMpegWriter

import numpy as np
from numpy.typing import NDArray, ArrayLike
from dataclasses import dataclass
from typing import override

from .util import split
from .io import init_input_filestream, close_filestream
from .constants import DATA_FILE_NAME, ANIMATION_FILE_NAME, \
    RADIAL_DENSITY_DISTRIBUTION_FILE_NAME


class ProcessingModule(ABC):
    r"""
    This is an abstract class for all "processing modules" so that they can
    operate in parallel to reduce the amount of data read to the disk at one
    time (such as generating figures, performing calculations, etc).
    """
    output_dir: Path

    def __init__(self, output_dir: Path):
        """
        :param output_dir: The output directory to write to.
        :type  output_dir: Path
        """
        self.output_dir = output_dir

    def begin(self) -> None:
        """
        Instructs the processing module to initialize such that it can
        produce data.
        """
        pass

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
    particle_radius: int = 20
    orientation_width: int = 2
    orientation_length: int = 3

    # color scheme
    head_color: str = 'tab:red'
    tail_color: str = 'tab:blue'
    particle_color: str = 'tab:cyan'
    orientation_color: str = 'tab:orange'

    @override
    def __init__(self, output_dir: Path):
        super().__init__(output_dir)

    @override
    def begin(self,
              fps: int = 20,
              max_bounds: ArrayLike = [ 100, 100, 100 ],
              show_debug_text: bool = True,
              show_heads_and_tails: bool = False) -> None:
        r"""
        Sets up the MatPlotLib animation for rendering.
        """
        self.show_debug_text = show_debug_text
        self.show_heads_and_tails = show_heads_and_tails

        self.fig = plt.figure(figsize=(8, 8))
        self.ax = self.fig.add_subplot(projection='3d')
        self.frames = []

        self.padding = 0.1
        self.max_bounds = np.asarray(max_bounds, dtype=float)

        self.output_filename = self.output_dir / ANIMATION_FILE_NAME
        self.writer = FFMpegWriter(fps)
        self.writer.setup(self.fig, str(self.output_filename), dpi=100)

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
    def append_state(self, system: NDArray, time: float) -> None:
        r"""
        Draws the current state to the animation.
        """
        positions, orientations = split(system)

        self.ax.clear()
        self.set_limits()

        self.ax.scatter(
            positions[:, 0],
            positions[:, 1],
            positions[:, 2],
            s = self.particle_radius,
            c = self.particle_color
        )

        self.ax.quiver(
            positions[:, 0],
            positions[:, 1],
            positions[:, 2],
            orientations[:, 0],
            orientations[:, 1],
            orientations[:, 2],
            length=self.orientation_length,
            normalize=True,
            linewidth=self.orientation_width,
            color=self.orientation_color,
        )

        self.ax.set_title(f't={time:.2f}')
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
    def __init__(self, output_dir: Path):
        super().__init__(output_dir)

    @override
    def begin(self, fps: int = 30) -> None:
        self.fig, self.ax = plt.subplots(figsize=(6, 4))

        self.bins = np.linspace(0, self.max_radius, self.n_bins + 1)
        self.bar_container = None

        self.output_filename = self.output_dir / \
            RADIAL_DENSITY_DISTRIBUTION_FILE_NAME
        self.writer = plt.matplotlib.animation.FFMpegWriter(fps=10)
        self.writer.setup(self.fig, str(self.output_filename), dpi=100)

        self.ax.set_xlabel("Radius")
        self.ax.set_ylabel("Count")
        self.ax.set_title("Radial Density Distribution")

    @override
    def append_state(self, system: NDArray, time: float) -> None:
        positions, _ = split(system)

        # center of mass
        com = np.mean(positions, axis=0)

        # radii
        radii = np.linalg.norm(positions - com, axis=1)

        counts, _ = np.histogram(radii, bins=self.bins)

        self.ax.clear()

        bin_centers = 0.5 * (self.bins[:-1] + self.bins[1:])

        self.ax.bar(bin_centers, counts, width=self.bins[1] - self.bins[0])

        self.ax.set_xlim(0, self.max_radius)
        self.ax.set_ylim(0, max(1, np.max(counts)))
        self.ax.set_title(f"Radial Density (t={time:.2f})")
        self.ax.set_xlabel("Radius")
        self.ax.set_ylabel("Count")

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


def process_data(
        output_dir: Path,
        generate_animation: bool=True,
        generate_density_animation: bool=True) -> None:
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
    io = init_input_filestream(datafile)

    modules: list[ProcessingModule] = []

    if generate_animation:
        modules.append(AnimationGenerator(output_dir))

    if generate_density_animation:
        modules.append(DensityAnimationGenerator(output_dir))

    for module in modules:
        module.begin()

    try:
        for i in range(len(io.time_dataset)):

            system = io.state_dataset[i]
            time = io.time_dataset[i]

            for module in modules:
                module.append_state(system, time)

    finally:

        for module in modules:
            module.end()

    close_filestream(io)
