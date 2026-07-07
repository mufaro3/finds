# implementation of barnes-hut and the FMM in 2-D
from sys import exit
from typing import Optional

import numpy as np
from numpy.typing import NDArray
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import solve_ivp
from tqdm import tqdm
from numba import njit, prange, optional, int32, float32
from numba.experimental import jitclass
from numba.types import Array

NDIMENSIONS=2
BOUNDS=np.array([20,20])
ORIENTATION_LENGTH=1

# NOTE TO SELF
#
# i have to manually define getters and setters for the array in order for
# this to actually work with Numba and SciPy simultaneously FUCK MY LIFE

# setters have to work on a single position


START = 0
END   = 1

@njit(inline='always')
def index_range(start: int, width: int):
    return (start, start + width)

@njit(inline='always')
def index_range_extend(previous_range: tuple[int, int], width: int):
    return index_range(previous_range[END], width)

@njit(inline='always')
def index_range_width(range_tuple: tuple[int, int]):
    return range_tuple[END] - range_tuple[START]

POSITION_RANGE    = index_range(0, NDIMENSIONS)
ORIENTATION_RANGE = index_range_extend(POSITION_RANGE, NDIMENSIONS)
LENGTH_RANGE      = index_range_extend(ORIENTATION_RANGE, 1)
VOL_FLOW_RANGE    = index_range_extend(LENGTH_RANGE, 1)
MASS_RANGE        = index_range_extend(VOL_FLOW_RANGE, 1)

# metadata

DATA_SLOTS = index_range_width(LENGTH_RANGE) + \
        index_range_width(VOL_FLOW_RANGE) + index_range_width(MASS_RANGE)
FISH_SLOTS = index_range_width(POSITION_RANGE) + \
        index_range_width(ORIENTATION_RANGE) + DATA_SLOTS

@njit(inline='always')
def fish_range_from_index(fish_index: int):
    return index_range(fish_index * FISH_SLOTS, FISH_SLOTS)

@njit(inline='always')
def subrange(fish_index: int, search_range: tuple[int, int]):
    return index_range(
        fish_index * FISH_SLOTS + search_range[START],
        fish_index * FISH_SLOTS + search_range[END]
    )

@nijt(inline='always')
def system_size(total_slots: int) -> int:
    return total_slots // FISH_SLOTS

# feature positions

FRONT_RANGE = index_range(0, NDIMENSIONS)
BACK_RANGE  = index_range_extend(FRONT_RANGE, NDIMENSIONS)
DATA_RANGE  = index_range_extend(BACK_RANGE, DATA_SLOTS)

assert index_range_width(FRONT_RANGE) + index_range_width(BACK_RANGE) + \
        DATA_SLOTS == FISH_SLOTS

@njit(inline='always')
def system_get_positions(system: NDArray) -> NDArray:
    return system[POSITION_RANGE[START]:POSITION_RANGE[END]]

@njit(inline='always')
def system_get_position_at(system: NDArray, index: int = None) -> NDArray:
    index_range = fish_range_from_index(index)
    return system[POSITION_RANGE[START]:POSITION_RANGE[END]]

@njit(inline='always')
def system_set_positions(system: NDArray, position_data: NDArray):
    system[:, POSITION_RANGE[START]:POSITION_RANGE[END]] = data

mocksystem = np.concatenate([
    [1, 2, 3, 4, 5, 6, 7, 8],
    [10, 20, 30, 40, 50, 60, 70, 80]
])
print(system_get_position_at(mocksystem, 1))
exit()

"""
def system_set_positions()
def system_get_orientations()
def system_set_orientations()
def system_get_lengths()
def system_set_lengths()
def system_get_volumetric_flow_rates()
def system_set_volumetric_flow_rates()
def system_get_masses()
def system_set_masss()
"""

def generate_system(N):
    system = np.empty((N,FISH_SLOTS), dtype=np.float32)
    system[:,0:NDIMENSIONS] = np.random.uniform(-BOUNDS, BOUNDS, size=(n,NDIMENSIONS))
    system[:,NDIMENSIONS:NDIMENSIONS*2] = -system.position / \
        np.linalg.norm(system.position, axis=1, keepdims=True)
    system.lengths[:] = np.random.uniform(0.5, 1.0, size=n)
    system.vol_flow_rates[:] = np.random.uniform(5, 20, size=n)
    system.masses[:] = np.random.uniform(0.1, 1.0, size=n)
    return system

def plot_system(system):
    plt.figure()

    endpoints = system.positions + ORIENTATION_LENGTH * system.orientations
    for p, e in zip(system.positions, endpoints):
        start = [p[0], e[0]]
        end = [p[1], e[1]]
        plt.plot(start, end, color='black')

    plt.scatter(*system.positions.T, color='red')

    ax = plt.gca()
    ax.set_xlim([-BOUNDS[0],BOUNDS[0]])
    ax.set_ylim([-BOUNDS[1],BOUNDS[1]])

    plt.xlabel('X')
    plt.ylabel('Y')

    plt.show()

@njit
def calculate_interaction(position_a, position_b):
    displacement = position_a - position_b
    return displacement / (np.linalg.norm(displacement) * NDIMENSIONS)

FEATURE_DTYPE = [
    ('size',  int32),
    ('front', float32[:, :]),
    ('back',  float32[:, :])
]

@jitclass(FEATURE_DTYPE)
class Features:
    def __init__(self, n):
        self.size  = n
        self.front = np.empty((n, NDIMENSIONS), dtype=np.float32)
        self.back  = np.empty((n, NDIMENSIONS), dtype=np.float32)
    def show(self):
        print(self.front)
        print(self.back)

@njit
def calculate_velocity_contribution(
        feature_positions, system, index_a, index_b):

    def interaction_to_velocity(interaction):
        return system[index_b, VOL_FLOW] / (4 * np.pi) * interaction

    external_velocities = Features(system.size)

    front_front = calculate_interaction(
        feature_positions.front[index_a],
        feature_positions.front[index_b]
    )

    front_back = calculate_interaction(
        feature_positions.front[index_a],
        feature_positions.back[index_b]
    )

    external_velocities.front = \
        interaction_to_velocity(front_front - front_back)

    back_front = calculate_interaction(
        feature_positions.back[index_a],
        feature_positions.front[index_b]
    )

    back_back = calculate_interaction(
        feature_positions.back[index_a],
        feature_positions.back[index_b]
    )

    external_velocities.back = \
        interaction_to_velocity(back_front - back_back)

    return external_velocities

@njit
def calculate_feature_positions(system):
    delta = system.length * system.orientation / 2
    feature_positions = Features(system.size)

    feature_positions.front = system.positions + delta
    feature_positions.back  = system.positions - delta

    return feature_positions

@njit(parallel=True)
def compute_external_velocity_brute_force(system):
    feature_positions = calculate_feature_positions(system)

    external_velocity_contribution = np.zeros((system.size, NDIMENSIONS*2))
    for i in prange(N):
        for j in range(N):
            if i == j:
                continue
            external_velocity_contribution[i] += \
                calculate_velocity_contribution(
                    feature_positions, system, i, j)

    return external_velocity_contribution

@njit
def compute_external_velocity_barnes_hut(system, theta):
    pass

@njit
def compute_external_velocity_fmm(system):
    pass

@njit
def calculate_feature_velocity(system, method, theta):
    internal_contribution = system.vol_flow_rates / \
        (4 * np.pi * system.lengths ** 2) * system.orientations
    internal_contribution = \
        np.hstack((internal_contribution, internal_contribution))

    external_contribution = None
    match method:
        case "brute":
            external_contribution = \
                compute_external_velocity_brute_force(system)
        case "barnes":
            external_contribution = \
                compute_external_velocity_barnes_hut(system, theta)
        case "FMM":
            external_contribution = \
                compute_external_velocity_fmm(system)
        case _:
            raise ValueError('method must be brute, barnes, or FMM')

    return internal_contribution + external_contribution

@njit
def calculate_derivative(system, method, theta):
    feature_velocity = calculate_feature_velocity(system, method, theta)
    velocity_diff = feature_velocity.front - feature_velocity.back

    # translational velocity
    translational_velocity = \
        (feature_velocity.front + feature_velocity.back) / 2

    # rotational velocity
    lagrange_multiplier = \
        -np.sum(velocity_diff * orientations, axis=1)[:, np.newaxis] / 2
    rotational_velocity = (velocity_diff + 2 * lagrange_multiplier *
                           orientations) / system.lengths

    return translational_velocity, rotational_velocity

def calculate_time_progression(system, method, theta,
                               rtol, atol, integration_method,
                               time_step, end_time, show_progress):
    progress_bar = tqdm(
        total      = end_time,
        disable    = not show_progress,
        desc       = 'Time Progression',
        bar_format = '{l_bar}{bar}| {n:.2f}/{total_fmt} s [{elapsed}]'
    )

    def differential_equation_functor(t, y):
        progress_bar.update(t - progress_bar.n)

        translational, rotational = calculate_derivative(
            y.reshape(system.shape), method, theta)

        return np.concatenate((translational.ravel(), rotational.ravel()))

    solved_state = solve_ivp(
        fun    = differential_equation_functor,
        t_span = (0, end_time),
        y0     = system.toarray(),
        method = integration_method,
        rtol   = rtol,
        atol   = atol,
        t_eval = np.arange(0, end_time, time_step)
    )

    times = solved_state.t
    trajectory = solved_state.y.T.reshape(-1, *system.shape)

    return trajectory, times

def animate_system(trajectory, times):
    fig, ax = plt.subplots()

    system0 = trajectory[0]
    positions0 = system0[:, POSITION]
    endpoints0 = positions0 + ORIENTATION_LENGTH * system0[:, ORIENTATION]

    lines = [
        ax.plot([p[0], e[0]], [p[1], e[1]], color='black')[0]
        for p, e in zip(positions0, endpoints0)
    ]

    scatter = ax.scatter(*positions0.T, color='red')

    ax.set_xlim([-BOUNDS[0], BOUNDS[0]])
    ax.set_ylim([-BOUNDS[1], BOUNDS[1]])
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect("equal")

    def update(frame):
        system = trajectory[frame]
        positions = system[:, POSITION]
        orientations = system[:, ORIENTATION]
        endpoints = positions + ORIENTATION_LENGTH * orientations

        ax.set_title(f"t={times[frame]:.2f}")
        for line, p, e in zip(lines, positions, endpoints):
            line.set_data(
                [p[0], e[0]],
                [p[1], e[1]]
            )

        scatter.set_offsets(positions)
        return [scatter, *lines]

    animation = FuncAnimation(
        fig, update, frames=times.size, interval=20, blit=True)

    plt.show()

result = calculate_time_progression(
    system = generate_system(10),
    method = 'brute',
    theta = 0,
    rtol = 1e-3,
    atol = 1e-6,
    integration_method = 'RK45',
    time_step = 1e-2,
    end_time = 10,
    show_progress = True
)

animate_system(*result)
