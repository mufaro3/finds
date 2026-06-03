import numpy as np
from matplotlib import pyplot as plt
import numba
from numpy.typing import NDArray, ArrayLike
from dataclasses import dataclass
from pprint import pprint
from types import SimpleNamespace
from numbers import Number

FISH_LENGTH=1 # m
VOLUMETRIC_FLOW_RATE=4*np.pi # rad/s
FISH_ORIENTATION_DELTA=0.01 # rad
NUMBER_OF_FISH=3
FISH_GENERATION_BOUNDS=[10,10,10]
FISH_SELF_PROPELLED_SPEED=VOLUMETRIC_FLOW_RATE/(4*np.pi*FISH_LENGTH**2)

TIME_STEP=0.01
END_TIME=1

OUTPUT_FILENAME='output/fish_movement.csv'

@dataclass
class SystemState:
    positions: NDArray[float]
    orientations: NDArray[float]

    def __add__(self, other):
        if isinstance(other, SystemState):
            return SystemState(
                self.positions + other.positions,
                self.orientations + other.orientations
            )
        else:
            raise ValueError(
                f'Only addition with SystemState is allowed. Given type {type(other)}')

    def __sub__(self, other):
        if isinstance(other, SystemState):
            return SystemState(
                self.positions + other.positions,
                self.orientations + other.orientations
            )
        else:
            raise ValueError(
                f'Only subtraction with SystemState is allowed. Given type {type(other)}')

    def __mul__(self, other):
        if isinstance(other, Number):
            return SystemState(
                self.positions * other,
                self.orientations * other
            )
        else:
            raise ValueError(
                f'Only division with a scalar is allowed. Given type {type(other)}')

    def __div__(self, other):
        if isinstance(other, Number):
            return SystemState(
                self.positions / other,
                self.orientations / other
            )
        else:
            raise ValueError(
                f'Only division with a scalar is allowed. Given type {type(other)}')
        
def spherical_to_cartesian(vec_spherical: NDArray) -> NDArray:
    """Converts a 2D spherical unit vector to a 3D cartesian vector"""
    theta, phi = vec_spherical[..., 0], vec_spherical[..., 1]
    cartesian = \
        np.stack([
            np.sin(phi) * np.cos(theta),
            np.sin(phi) * np.sin(theta),
            np.cos(phi)
        ], axis=0)
    return cartesian
    
def generate_fish(bounds: ArrayLike, angle_delta: float) -> NDArray:
    if angle_delta >= np.pi:
        raise ValueError(
            f'Angular perturbation is too large! Maximum of pi radians. angle_delta={angle_delta::2e} rad.')
    
    bounds = np.asarray(bounds)
    position = np.random.uniform(-bounds, bounds, bounds.shape)
    orientation = np.random.uniform(-angle_delta, angle_delta, 2)
    # NOTE: orientations is defined as [theta, phi]
    
    return position, spherical_to_cartesian(orientation)

def generate_school(n: int, bounds: ArrayLike, angle_delta: float = 0) -> NDArray:
    """Generates a fresh school/state at random"""
    raw = [generate_fish(bounds, angle_delta) for i in range(n)]
    positions, orientations = map(np.stack, zip(*raw))
    return SystemState(positions, orientations)

def calculate_feature_positions(school: SystemState) -> tuple[NDArray]:
    """Computes the head and tail positions for the school"""
    delta = school.orientations * FISH_LENGTH / 2
    heads = school.positions + delta
    tails = school.positions - delta
    return heads, tails

def compute_pairwise_interactions(
        school: SystemState,
        feature_positions: tuple[NDArray],
        use_barnes_hut: bool,
        head: bool) -> NDArray:
    """Computes the sum of the pairwise interactions for each fish"""
    pass

def calculate_feature_velocities(
        school: SystemState,
        feature_positions: tuple[NDArray],
        use_barnes_hut: bool) -> NDArray:
    """Calculates the velocities for the head and tail of all fish in a system"""
    internal_contrib = FISH_SELF_PROPELLED_SPEED * orientations
    pairwise_interactions_sum = np.array([
        compute_pairwise_interactions(
            school, feature_positions, use_barnes_hut, head=True),
        compute_pairwise_interactions(
            school, feature_positions, use_barnes_hut, head=False)
    ])
    external_contrib = VOLUMETRIC_FLOW_RATE / (4 * np.pi) *\
        pairwise_interactions_sums
    return internal_contrib + external_contrib
    
def calculate_state_derivative(
        school: SystemState, use_barnes_hut: bool) -> SystemState:
    """Computes the derivative of the system (school) matrix"""
    feature_positions = calculate_feature_positions(school)

    head_velocities, tail_velocities = calculate_feature_velocities(
        school,
        feature_positions,
        use_barnes_hut,
        head=False
    )

    translational_derivative = \
        (head_velocities + tail_velocities) / 2

    velocity_diff = head_velocities - tail_velocities
    # computing the dot product as a matrix multiplication for speed
    lagrange_mult = -np.einsum('ij,ij->i', velocity_diff, school.orientations) / 2
    rotational_derivative = \
        (velocity_diff + 2 * lagrange_mult * school.orientations) / FISH_LENGTH

    return SystemState(
        translational_derivative,
        rotational_derivative
    )
    
def calculate_update_rk4(school: SystemState, time_step: float) -> SystemState:
    """Computes the updated state of the school after the next time step using RK4"""
    # some simple aliases to make the calculations more mathematical
    X0 = school
    f  = calculate_state_derivative
    dt = time_step

    k1 = f(X0)
    k2 = f(X0 + k1 * dt / 2)
    k3 = f(X0 + k2 * dt / 2)
    k4 = f(X0 + k3 * dt)

    k = (k1 + 2 * k2 + 2 * k3 + k4) / 6

    Xt = X0 + k * dt

    updated_state = Xt

    # normalizes the orientation vectors
    updated_state.orientations = updated_state.orientations /\
        np.linalg.norm(updated_state.orientations, axis=1, keepdims=True)
    
    return updated_state
    
def serialize_to_file(school: SystemState, time: float, output_filename: str) -> None:
    """Writes the system state to the datafile"""
    pass

def perform_simulation():
    school = generate_school(
        NUMBER_OF_FISH,
        FISH_GENERATION_BOUNDS,
        FISH_ORIENTATION_DELTA
    )
    pprint(school)

    # serialize the starting state
    serialize_to_file(school, time, OUTPUT_FILENAME)
    
    return
    
    time=0
    while time<END_TIME:
        school = calculate_update_rk4(school, TIME_STEP)
        serialize_to_file(school, time, OUTPUT_FILENAME)
        time += TIME_STEP
    
if __name__ == '__main__':
    perform_simulation()
