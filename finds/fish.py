import numpy as np
from numba import njit
from numpy.typing import ArrayLike, NDArray

from .util import rejoin, split


@njit
def normalize_orientation_vectors(system: NDArray) -> NDArray:
    r"""
    Normalizes the orientation vectors for the given system.

    :param system: The system to normalize.
    :type system:  NDArray

    :returns: The system, with normalized orientation vectors.
    :rtype:   NDArray
    """
    pos, ori = split(system)
    norms = np.sqrt(np.sum(ori**2, axis=1))

    normalized_ori = ori / norms[:, np.newaxis]
    return rejoin(pos, normalized_ori)


def perturb_orientations(orientations: NDArray, angle_delta: float) -> NDArray:
    r"""
    Perturbs an orientations matrix by :math:`\Delta \theta`.

    :param orientations: The matrix of fish orientations
    :type  orientations: NDArray

    :param angle_delta: The maximum random angular perturbation for a given
      fish in radians, :math:`\Delta \theta`
    :type angle_delta: float

    :rtype: NDArray
    """
    if orientations is None:
        raise ValueError('Orientations is NoneType.')

    if orientations.ndim != 2 or orientations.shape[1] != 3:
        raise ValueError('Orientations matrix should have shape (N,3).')

    if angle_delta >= np.pi or angle_delta < 0:
        raise ValueError("Angular perturbation is too large! "+
                         "angle_delta=" + str(angle_delta))

    for i in range(orientations.shape[0]):
        theta = np.random.uniform(-angle_delta, angle_delta)
        phi = np.random.uniform(-angle_delta, angle_delta)

        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)

        sin_phi = np.sin(phi)
        cos_phi = np.cos(phi)

        Rz = np.array([
            [1, 0,          0        ],
            [0, cos_theta, -sin_theta],
            [0, sin_theta,  cos_theta]
        ])

        Rx = np.array([
            [cos_phi, -sin_phi, 0],
            [sin_phi,  cos_phi, 0],
            [0,        0,       1]
        ])

        orientations[i] = Rz @ Rx @ orientations[i]

    return orientations


def generate_positions_random(n: int, bounds: ArrayLike) -> NDArray:
    r"""
    Generates a fresh system at random.

    :param n: The number of fish to include in the system.
    :type n: int

    :param bounds: The cartesian bounds (in meters) for which to generate fish.
      Supplied in form :math:`(x_b,y_b,z_b)`, and means that the fish will be
      generated at position :math:`(x,y,z)` where :math:`x \in (-x_b, x_b),
      y \in (-x_b, y_b)`, and :math:`z \in (-z_b, z_b)`.
    :type bounds: ArrayLike

    :param angle_delta: The maximum random angular perturbation for all fish
      in radians.
    :type angle_delta: float

    :returns: A matrix with shape :math:`(N,6)` storing the position of
      each fish.
    :rtype: NDArray
    """
    bounds = np.asarray(bounds)
    if bounds.shape != (3,):
        raise ValueError('Bounds must be 3 dimensional.')

    positions = np.random.uniform(-bounds, bounds, size=(n, 3))
    return positions


def generate_positions_lattice(side_length: int, spacing: float) -> NDArray:
    r"""
    Generates the positions of fish arranged as a cubic lattice centered
    at (0,0).

    :param side_length: The side length of the cube.
    :param spacing: The grid spacing between fish within the lattice.
    :rtype: NDArray
    """
    half = side_length / 2

    x = np.arange(-half, half + spacing, spacing)
    y = np.arange(-half, half + spacing, spacing)
    z = np.arange(-half, half + spacing, spacing)

    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

    positions = np.column_stack([
        X.ravel(),
        Y.ravel(),
        Z.ravel()
    ])

    return positions


def generate_positions_sphere(radius: int, spacing: float) -> NDArray:
    r"""
    Generates the positions of fish arranged as a sphere centered at (0,0).

    :param radius:  The radius of the cube.
    :param spacing: The grid spacing between fish within the sphere.
    :rtype NDArray:
    """
    r = radius

    x = np.arange(-r, r + spacing, spacing)
    y = np.arange(-r, r + spacing, spacing)
    z = np.arange(-r, r + spacing, spacing)

    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

    R2 = X**2 + Y**2 + Z**2
    mask = R2 <= r**2

    positions = np.column_stack([
        X[mask],
        Y[mask],
        Z[mask]
    ])

    return positions


def generate_positions_square(side_length: int, spacing: float) -> NDArray:
    r"""
    Generates the positions of fish arranged as a planar square oriented with
    its normal vector on the :math:`x`-axis and centered at (0,0).

    :param side_length: The side length of the square.
    :param spacing: The grid spacing between fish in the square.
    :rtype: NDArray
    """
    half = side_length / 2

    y = np.arange(-half, half + spacing, spacing)
    z = np.arange(-half, half + spacing, spacing)

    Y, Z = np.meshgrid(y, z, indexing="ij")

    x = np.zeros_like(Y)

    positions = np.column_stack([
        x.ravel(),
        Y.ravel(),
        Z.ravel()
    ])

    return positions


def generate_positions_circle(radius: int, spacing: float):
    r"""
    Generates the positions of fish arranged as a planar circle oriented with
    its normal vector on the :math:`x`-axis and centered at (0,0).

    :param radius: The radius of the circle
    :param spacing: The grid spacing between fish in the circle.
    :rtype: NDArray
    """
    r = radius

    y = np.arange(-r, r + spacing, spacing)
    z = np.arange(-r, r + spacing, spacing)

    Y, Z = np.meshgrid(y, z, indexing="ij")

    mask = (Y**2 + Z**2) <= r**2

    positions = np.column_stack([
        np.zeros(np.sum(mask)),
        Y[mask],
        Z[mask]
    ])

    return positions


def generate_system(
        distribution: str,
        orientation: str,
        size: int = 20,
        spacing: float = 1.0,
        n_random: int = 2,
        bounds: ArrayLike = [10, 10, 10],
        angle_delta: float = 0) -> NDArray:
    r"""
    Generates a system of fish distributed according to :code:`distribution`:
    either spherically, as a lattice, or randomly within a box.

    :param distribution: How the fish are distributed.
    :type  distribution: 'lattice' | 'sphere' | 'random' | 'square' | 'circle'

    :param orientation: How the fish are oriented.
    :type  orientation: 'aligned' | 'radial outward' | 'radial inward'

    :param n_random: The number of fish, provided that this is a random
      distribution. For all other distributions, the :code:`size` parameter is
      used.

    :type  bounds: ArrayLike

    :param spacing: How far apart the fish are (assuming it isn't random).
    :type  spacing: float

    :rtype: NDArray
    """
    positions = None
    orientations = None

    if distribution == 'random':
        positions = generate_positions_random(
            n=n_random,
            bounds=bounds
        )
    elif distribution == 'lattice':
        positions = generate_positions_lattice(
            side_length=size,
            spacing=spacing
        )
    elif distribution == 'sphere':
        positions = generate_positions_sphere(
            radius=size,
            spacing=spacing
        )
    elif distribution == 'square':
        positions = generate_positions_square(
            side_length=size,
            spacing=spacing
        )
    elif distribution == 'circle':
        positions = generate_positions_circle(
            radius=size,
            spacing=spacing
        )
    else:
        raise ValueError(
            f'Invalid distribution parameter \"{distribution}\". '+\
            'Must be one of: lattice, sphere, random, square, circle')

    N = positions.shape[0]

    if orientation == 'aligned':
        orientations = \
            np.tile(np.array([1.0, 0.0, 0.0], dtype=np.float64), (N, 1))
    elif orientation == 'radial outward':
        orientations = positions.copy()
    elif orientation == 'radial inward':
        orientations = -positions.copy()
    else:
        raise ValueError(f'Invalid orientation parameter \"{orientation}\". '+
                         'Must be one of: \"aligned\", \"radial outward\", '+
                         '\"radial inward\".')

    if angle_delta != 0:
        orientations = perturb_orientations(orientations, angle_delta)

    # remove any zeros
    norms = np.linalg.norm(orientations, axis=1, keepdims=True)
    zero_mask = norms.squeeze() == 0
    orientations[zero_mask] = np.array([1.0, 0.0, 0.0])

    system_raw = rejoin(positions, orientations)

    return normalize_orientation_vectors(system_raw)
