import numpy as np
from numba import njit
from numpy.typing import NDArray


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
    if mat is None:
        raise ValueError('Matrix is NoneType.')

    if mat.ndim == 1:
        vec = mat

        if vec.size != 6:
            raise ValueError('Vector does not contain 6 elements.')

        return vec[:3], vec[3:]

    if mat.ndim != 2:
        raise ValueError('Expected 2D matrix.')
    if mat.shape[1] != 6:
        raise ValueError("Expected matrix with shape (N,6).")

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
    return np.concatenate((tophalf, bottomhalf), axis=-1)
