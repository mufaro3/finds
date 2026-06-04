Non-Mathematical Fish/System Management Functions
=================================================

For this simulation, the system of fish are stored the following format:

.. math::

   \mathbf{X} := \begin{bmatrix}
   x_1 & y_1 & z_1 & n_{1x} & n_{1y} & n_{1z} \\
   x_2 & y_2 & z_2 & n_{2x} & n_{2y} & n_{2z} \\
   \vdots & \vdots & \vdots & \vdots & \vdots & \vdots \\
   x_N & y_N & z_N & n_{Nx} & n_{Ny} & n_{Nz}
   \end{bmatrix}

This module consists of fish management functions that aren't specified in papers
:cite:t:`mabrouk2024` and :cite:t:`mabrouk2025`, such as :py:func:`src.fish.positions`, :py:func:`src.fish.orientations`, and :py:func:`src.fish.normalize_orientation_vectors` alongside fish/system generation functions.

.. automodule:: src.fish
   :members:
   :undoc-members:
   :show-inheritance:
