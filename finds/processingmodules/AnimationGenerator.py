from pathlib import Path
from typing import override

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FFMpegWriter
from numpy.typing import ArrayLike, NDArray

from .ProcessingModule import ProcessingModule

from ..util import split


class AnimationGenerator(ProcessingModule):
    r"""
    Generates a 3-D animation of the full fish simulation, taking (0,0) to be
    the origin. Allows for rapid customization of simulation colors, viewing,
    and debug parameters as well.

    Attributes
    ==========
    fps: int
      The frames-per-second of the simulation, or how fast the animation will
      be. Default is 30 FPS.

    dpi: int
      The resolution of the video in dots per inch. Defaults to 100.

    video_output_filename: str

    :param max_bounds: bruh
    :type max_bounds: ArrayLike
    """
    def __init__(
            self, *,
            # video quality and information
            fps: int = 30,
            dpi: int = 100,
            video_output_filename: str = 'animation3d.mp4',

            # sizing
            max_bounds: ArrayLike = [ 100, 100, 100 ],
            particle_radius: int = 30,
            orientation_width: int = 5,
            orientation_length: int = 5,
            padding: float = 0.1,
            figsize: tuple[int] = (8,8),

            # color scheme
            head_color: str = 'tab:red',
            tail_color: str = 'tab:blue',
            particle_color: str = 'tab:cyan',
            orientation_color: str = 'tab:orange'
    ):
        self.fps = fps
        self.dpi = dpi
        self.video_output_filename = video_output_filename

        self.max_bounds = np.asarray(max_bounds, dtype=float)
        self.particle_radius = particle_radius
        self.orientation_width = orientation_width
        self.orientation_length = orientation_length
        self.padding = padding
        self.figsize = figsize

        self.head_color = head_color
        self.tail_color = tail_color
        self.particle_color = particle_color
        self.orientation_color = orientation_color

    @override
    def begin(self, output_dir: Path) -> None:
        r"""
        Sets up the MatPlotLib animation for rendering.

        :type max_bounds: ArrayLike
        """
        super().begin(output_dir)

        # set the figure and axis
        self.fig = plt.figure(figsize=self.figsize)
        self.ax = self.fig.add_subplot(projection='3d')

        # set the limits
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

        self.ax.view_init(elev=20, azim=45)

        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        self.ax.set_zlabel("Z")

        # initialize the position particles
        self.scatter_plot = self.ax.scatter(
            [], [], [],
            s=self.particle_radius,
            c=self.particle_color
        )

        # initialize the orientations
        self.quiver_plot = self.ax.quiver(
            0, 0, 0, 0, 0, 0,
            color=self.orientation_color
        )

        # define information for the video (title, writer)
        self.title_text = self.ax.set_title('')
        self.output_filename = self.output_dir / self.video_output_filename
        self.writer = FFMpegWriter(self.fps)
        self.writer.setup(self.fig, str(self.output_filename), dpi=self.dpi)

    @override
    def append_state(self, system: NDArray, time: float) -> None:
        r"""
        Draws the current state to the animation.
        """
        positions, orientations = split(system)

        # scatterplot the positions in 3D
        self.scatter_plot._offsets3d = (
            positions[:, 0],
            positions[:, 1],
            positions[:, 2]
        )

        # update the orientation vectors
        vecs = orientations * self.orientation_length
        self.quiver_plot.set_segments(
            np.stack([positions, positions + vecs], axis=1)
        )

        # set the title
        self.title_text.set_text(f't={time:.2f}')

        # update the frame
        self.writer.grab_frame()

    @override
    def end(self) -> None:
        r"""
        Closes the animation and saves the file.
        """
        self.writer.finish()
        plt.close(self.fig)
        print(f'Animation saved at {self.output_filename}')
