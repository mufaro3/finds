from numba import jit, njit
import numpy as np
from numpy.typing import NDArray, ArrayLike
from .util import spherical_to_cartesian
    
@njit
def generate_fish(bounds: ArrayLike, angle_delta: float) -> NDArray:
    r"""
    Generates a fish with a random position with :code:`bounds` and a random
    angular perturbation as specified by :code:`angle_delta`.

    :param bounds: The cartesian bounds (in meters) for which to generate fish.
      Supplied in form :math:`(x_b,y_b,z_b)`, and means that the fish will be
      generated at position :math:`(x,y,z)` where :math:`x \in (-x_b, x_b), y
      \in (-x_b, y_b)`, and :math:`z \in (-z_b, z_b)`.
    :type bounds: ArrayLike

    :param angle_delta: The maximum random angular perturbation for a given
      fish in radians.
    :type angle_delta: float

    :returns: An array in form :math:`[x, y, z, n_x, n_y, n_z]` storing the
      position and orientation of the fish.
    :rtype: NDArray
    """
    if angle_delta >= np.pi:
        raise ValueError("Angular perturbation is too large! "+
                         "angle_delta=" + str(angle_delta))
    
    bounds = np.asarray(bounds)

    position = np.zeros(bounds.shape[0], dtype=np.float64)
    for i in range(bounds.shape[0]):
        position[i] = np.random.uniform(-bounds[i], bounds[i])
        
    orientation = np.zeros(2, dtype=np.float64)
    for i in range(2):
        orientation[i] = np.random.uniform(0, angle_delta)
    orientation_cartes = spherical_to_cartesian(orientation)
        
    new_fish = np.zeros(6, dtype=np.float64)
    new_fish[:3] = position
    new_fish[3:] = orientation_cartes
    
    return new_fish

@njit
def generate_system(n: int, bounds: ArrayLike, angle_delta: float = 0) -> NDArray[float]:
    r"""
    Generates a fresh system at random.

    :param n: The number of fish to include in the system.
    :type n: int

    :param bounds: The cartesian bounds (in meters) for which to generate fish.
      Supplied in form :math:`(x_b,y_b,z_b)`, and means that the fish will be
      generated at position :math:`(x,y,z)` where :math:`x \in (-x_b, x_b),
      y \in (-x_b, y_b)`, and :math:`z \in (-z_b, z_b)`.
    :type bounds: ArrayLike

    :param angle_delta: The maximum random angular perturbation for all fish in radians.
    :type angle_delta: float

    :returns: A matrix with shape (N,6) storing the position and orientation of each fish.
    :rtype: NDArray
    """
    system = np.empty((n, 6), dtype=np.float64)

    for i in range(n):
        system[i] = generate_fish(bounds, angle_delta)

    return system

@njit
def normalize_orientation_vectors(system: NDArray) -> NDArray:
    r"""
    Normalizes the orientation vectors for the given system.

    :param system: The system to normalize.
    :type system:  NDArray

    :returns: The system, with normalized orientation vectors.
    :rtype:   NDArray
    """
    pos   = positions(system)
    ori   = orientations(system)
    norms = np.sqrt(np.sum(ori**2, axis=1))

    normalized = ori / norms[:, np.newaxis]
    recombined = np.hstack((pos, normalized))
    
    return recombined
