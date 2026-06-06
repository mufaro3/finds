from numba import jit, njit
import numpy as np
from numpy.typing import NDArray, ArrayLike

@njit
def split(mat: NDArray) -> tuple[NDArray, NDArray]:
    """
    Splits a matrix of shape (N,6) into two matrices of (N,3) and
    (N,3), or the top and bottom halves of the original matrix.

    :param mat: The original matrix of shape :math:`(N,6)`.
    :type  mat: NDArray

    :returns: The top and bottom halves of the matrix.
    :rtype: tuple[NDArray, NDArray]
    """
    if len(mat.shape) == 1:
        vec = mat
        if vec.size != 6:
            raise ValueError('Vector does not contain 6 elements.')
        return vec[:3], vec[3:]

    if len(mat.shape) != 2 or mat.shape[1] != 6:
        raise ValueError('Matrix must have shape (N,6).')

    return mat[:, :3], mat[:, 3:]

@njit
def rejoin(tophalf: NDArray, bottomhalf: NDArray) -> NDArray:
    """
    Rejoins the top and bottom halves of a system matrix of similar.
    Requires that :code:`tophalf` and :code:`bottomhalf` have the
    same shape, :math:`(N,6)`.

    :param tophalf: The top half of the matrix.
    :type  tophalf: NDArray

    :param bottomhalf: The bottom half of the matrix.
    :type  bottomhalf: NDArray

    :returns: The rejoined matrix.
    :rtype: NDArray
    """
    return np.hstack((tophalf, bottomhalf))

@njit
def iterate_excluding_self(array: NDArray):
    """
    Iterates over the array, excluding the self.

    :type array: NDArray
    :rtype: NDArray
    """
    for i, item in enumerate(array):
        others = np.concatenate((array[:i], array[i+1:]))
        yield item, others

@njit
def cartesian_to_spherical(vec_cartesian: ArrayLike) -> NDArray:
    r"""
    Converts  a 3D cartesian unit vector to a 2D spherical vector
    using the spherical to cartesian formulas

    .. math::
        :nowrap:
    
        \begin{align}
          \theta &= \arctan(y/x) \\
          \phi &= \arctan\left(\frac{\sqrt{x^2 + y^2}}{z}\right)
        \end{align}

    :param vec_cartesian: The cartesian vector.
    :type  vec_cartesian: ArrayLike
    
    :returns: The Cartesian vector in Spherical coordinates.
    :rtype: NDArray
    """
    vec_cartesian = np.asarray(vec_cartesian)

    if vec_cartesian.shape != (3,):
        raise ValueError('Cartesian vector is not three-dimensional!')
    
    norm = np.linalg.norm(vec_cartesian)
    if not np.isclose(norm, 1):
        raise ValueError('Vector ' + str(vec_cartesian) + ' is not a unit vector, '+
                         'norm = ' + str(norm))

    x = vec_cartesian[0]
    y = vec_cartesian[1]
    z = vec_cartesian[2]

    recovered_theta = np.atan2(y, x)
    recovered_phi   = np.atan2( np.sqrt( x ** 2 + y ** 2 ), z )

    # normalize back into bounds
    normalized_theta = np.mod(recovered_theta, 2 * np.pi)
    normalized_phi   = np.mod(recovered_phi,   np.pi)

    return np.array([normalized_theta, normalized_phi])
    
        
@njit
def spherical_to_cartesian(vec_spherical: ArrayLike) -> NDArray:
    r"""
    Converts a 2D spherical unit vector to a 3D cartesian vector
    using the spherical conversion formulas

    .. math::
        :nowrap:
    
        \begin{align}
        x &= \sin \phi \cos \theta \\
        y &= \sin \phi \sin \theta \\
        z &= \cos \phi
        \end{align}

    :param vec_spherical: A vector of spherical coordinates
      :math:`((\theta_1, \phi_1), (\theta_2, \phi_2), \dots, (\theta_n, \phi_n))`
    :type  vec_spherical: ArrayLike

    :returns: The vector in cartesian coordinates.
    :rtype: NDArray
    """
    vec_spherical = np.asarray(vec_spherical)
    
    if vec_spherical.shape != (2,):
        raise ValueError('Spherical Vector must have shape (2,), '+
                         'but was given shape: ' + str(vec_spherical.shape))
    
    theta = vec_spherical[0]
    phi   = vec_spherical[1]
    
    # verify that the input is actually angular
    if theta < 0 or theta >= 2 * np.pi:
        raise ValueError(f'Non-spherical argument. Theta={theta} is not within [0,2pi).')
    if phi < 0 or phi > np.pi:
        raise ValueError(f'Non-spherical argument. Phi={phi} is not within [0,pi].')
    
    cartesian = np.empty(3, dtype=np.float64)

    sin_phi = np.sin(phi)

    cartesian[0] = sin_phi * np.cos(theta)
    cartesian[1] = sin_phi * np.sin(theta)
    cartesian[2] = np.cos(phi)

    return cartesian
