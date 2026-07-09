# implementation of barnes-hut and the FMM in 2-D
from sys import exit
from typing import Optional

import numpy as np
from numpy.typing import NDArray
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import solve_ivp
from tqdm import tqdm
from numba import njit, prange, optional, int64, float64, float64, boolean, \
     uint64
from numba.experimental import jitclass
from numba.types import Array

NDIM=2
BOUNDS=np.array([200,200])
ORIENTATION_LENGTH=1

LENGTH_RANGE=(5,20)
VOLUMETRIC_FLOW_RANGE=(5,30)

# fish type
# position (3) orientation (3) length volumetric-flow

POSITION=slice(0,NDIM)
ORIENTATION=slice(NDIM,NDIM*2)
LENGTH=slice(NDIM*2,NDIM*2+1)
VOL_FLOW_RATE=slice(NDIM*2+1,NDIM*2+2)

DATA_SLOTS=2
FISH_SLOTS=2*NDIM+DATA_SLOTS

@njit
def norm(arr):
    if arr.ndim == 1:
        return np.linalg.norm(arr)

    if arr.ndim == 2:
        norms = np.empty(arr.shape[0], dtype=arr.dtype)
        for i in range(arr.shape[0]):
            norms[i] = np.linalg.norm(arr[i, :])
        return norms

    raise ValueError("arr should be 1D or 2D")

@njit
def generate_system(N: int, bounds: NDArray = BOUNDS) -> NDArray:
    system = np.empty((N,FISH_SLOTS), dtype=np.float64)

    # positions
    for i in range(N):
        for j in range(NDIM):
            system[i,j] = np.random.uniform(-bounds[j], bounds[j])

    # orientations
    system[:,ORIENTATION] = -system[:,POSITION] / \
        norm(system[:,POSITION]).reshape(N,1)

    # length
    system[:,LENGTH] = np.random.uniform(
        LENGTH_RANGE[0],
        LENGTH_RANGE[1],
        size=(N,1)
    )

    # volumetric flow rate
    system[:,VOL_FLOW_RATE] = np.random.uniform(
        VOLUMETRIC_FLOW_RANGE[0],
        VOLUMETRIC_FLOW_RANGE[1],
        size=(N,1)
    )

    assert np.isfinite(system).all() and np.isreal(system).all()

    return system

def plot_system(system):
    plt.figure()

    endpoints = system[:,POSITION] + ORIENTATION_LENGTH * system[:,ORIENTATION]
    for p, e in zip(system[:,POSITION], endpoints):
        start = [p[0], e[0]]
        end = [p[1], e[1]]
        plt.plot(start, end, color='black')

    plt.scatter(*system[:,POSITION].T, color='red')

    ax = plt.gca()
    ax.set_xlim([-BOUNDS[0],BOUNDS[0]])
    ax.set_ylim([-BOUNDS[1],BOUNDS[1]])

    plt.xlabel('X')
    plt.ylabel('Y')

    plt.show()

@njit
def calculate_interaction(position_a, position_b):
    displacement = position_a - position_b
    return displacement / (norm(displacement) ** NDIM)

# features matrix
FRONT = slice(0,NDIM)
BACK  = slice(NDIM,NDIM*2)

FEATURE_SLOTS = 4

@njit
def calculate_velocity_contribution(
        feature_positions, system, index_a, index_b):

    def interaction_to_velocity(interaction):
        return system[index_b, VOL_FLOW_RATE] / (4 * np.pi) * interaction

    external_velocities = np.empty((system.shape[0],FEATURE_SLOTS))

    front_front = calculate_interaction(
        feature_positions[index_a, FRONT],
        feature_positions[index_b, FRONT]
    )

    front_back = calculate_interaction(
        feature_positions[index_a, FRONT],
        feature_positions[index_b, BACK]
    )

    external_velocities[:,FRONT] = \
        interaction_to_velocity(front_front - front_back)

    back_front = calculate_interaction(
        feature_positions[index_a, BACK],
        feature_positions[index_b, FRONT]
    )

    back_back = calculate_interaction(
        feature_positions[index_a, BACK],
        feature_positions[index_b, BACK]
    )

    external_velocities[:,BACK] = \
        interaction_to_velocity(back_front - back_back)

    return external_velocities

@njit
def calculate_feature_positions(system):
    delta = system[:,LENGTH] * system[:,ORIENTATION] / 2
    feature_positions = np.empty((system.shape[0], FEATURE_SLOTS))

    feature_positions[:,FRONT] = system[:,POSITION] + delta
    feature_positions[:,BACK]  = system[:,POSITION] - delta

    return feature_positions

@njit(parallel=True)
def compute_external_velocity_brute_force(system):
    N = system.shape[0]
    feature_positions = calculate_feature_positions(system)
    external_velocity_contribution = np.zeros((N,FEATURE_SLOTS))

    for i in prange(N):
        for j in range(N):
            if i == j:
                continue
            external_velocity_contribution[i] += \
                calculate_velocity_contribution(
                    feature_positions, system, i, j)

    return external_velocity_contribution

@njit
def extend_1d(array, newshape):
    pass

@njit
def extend_2d(array, newshape):
    pass

@jitclass([
    ("center",           float64[:, :]),
    ("side_length",      float64[:]),

    ("morton",           uint64[:]),
    ("is_leaf",          boolean[:]),
    ("is_occupied",      boolean[:]),

    ("cluster",          float64[:, :]),
    ("cluster_size",     float64[:]),
    ("data",             float64[:, :])
    ("front_pos",        float64[:, :]),
    ("back_pos",         float64[:, :]),

    ("max_node_count", int64),
    ("total_stored_nodes", int64)
])
class Octree:
    def __init__(self, max_nodes):
        self.center      = np.empty((max_nodes, NDIM), dtype=np.float64)
        self.side_length = np.empty(max_nodes, dtype=np.float64)

        self.morton      = np.empty(max_nodes, dtype=np.float64)
        self.is_leaf     = np.empty(max_nodes, dtype=bool)
        self.is_occupied = np.zeros(max_nodes, dtype=bool)

        self.cluster      = np.empty((max_nodes, FISH_SLOTS), dtype=np.float64)
        self.cluster_size = np.empty(max_nodes, dtype=np.float64)
        self.data         = np.empty((max_nodes, FISH_SLOTS), dtype=np.float64)
        self.front_pos    = np.empty((max_nodes, NDIM), dtype=np.float64)
        self.back_pos     = np.empty((max_nodes, NDIM), dtype=np.float64)

        self.max_node_count = max_nodes
        self.total_stored_nodes = 0

    def extend(self):
        """double the amount of available node space"""
        self.max_node_count *= 2

        extend_2d(self.center,      shape=(self.max_node_count, NDIM))
        extend_1d(self.side_length, shape=self.max_node_count)

        extend_1d(self.morton,      shape=self.max_node_count)
        extend_1d(self.is_leaf,     shape=self.max_node_count)
        extend_1d(self.is_occupied, shape=self.max_node_count)

        extend_2d(self.cluster,      shape=(self.max_node_count, NDIM))
        extend_1d(self.cluster_size, shape=self.max_node_count)
        extend_2d(self.data,         shape=(self.max_node_count, NDIM))
        extend_2d(self.front_pos,    shape=(self.max_node_count, NDIM))
        extend_2d(self.back_pos,     shape=(self.max_node_count, NDIM))

    def sort(self):
        # sort according to the morton number
        pass

    def calculate_child_octant_index(self, parent_index, child_position):
        child_octant_index = 0
        for dim_idx in range(NDIM):
            if child_position[dim_idx] > self.center[parent_index, dim_idx]:
                child_octant_index += 2 ** dim_idx
        return child_octant_index

    def calculate_child_center_morton(self, parent_index, child_octant_index):
        offset_length = self.side_length[parent_index] / 4
        child_center = self.center[parent_index].copy()

        for dimension_index in range(NDIM):
            if child_octant_index & (2 ** dimension_bit):
                child_center[dimension_index] += offset_length
            else:
                center_center[dimension_index] -= offset_length

        child_morton = (self.morton[parent_index] <<< NDIM) | child_octant_index

    def insert(self, fish, node_index):
        if not self.is_occupied[node_index]:
            self.data[node_index] = fish
            self.is_occupied[node_index] = True
        elif self.is_leaf[node_index]

    def insert_into_children()

@njit
def build_octree(system, feature_positions):
    N = system.shape[0]
    octree = Octree(2 * N)
    for i in range(N):
        octree.insert(system[i], feature_positions[i])
    octree.sort()

@njit
def compute_external_velocity_singular_barnes_hut(octree, system, i):
    pass

@njit(parallel=True)
def compute_external_velocity_barnes_hut(system, theta):
    N = system.shape[0]
    feature_positions = calculate_feature_positions(system)
    external_velocity_contribution = np.zeros((N,FEATURE_SLOTS))
    octree = build_octree(system)

    for i in prange(N):
        external_velocity_contribution[i] = \
            compute_external_velocity_singular_barnes_hut(octree, system, i)

    return external_velocity_contribution

@njit
def compute_external_velocity_fmm(system):
    pass

@njit
def calculate_feature_velocity(system, method, theta):
    internal_contribution = system[:,VOL_FLOW_RATE] / \
        (4 * np.pi * system[:,LENGTH] ** 2) * system[:,ORIENTATION]
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
    velocity_diff = feature_velocity[:,FRONT] - feature_velocity[:,BACK]

    # translational velocity
    translational_velocity = \
        (feature_velocity[:,FRONT] + feature_velocity[:,BACK]) / 2

    # rotational velocity
    lagrange_multiplier = -np.sum(
        velocity_diff * system[:,ORIENTATION], axis=1
    )[:, np.newaxis] / 2
    rotational_velocity = (velocity_diff + 2 * lagrange_multiplier *
                           system[:,ORIENTATION]) / system[:,LENGTH]

    return np.hstack((
        translational_velocity,
        rotational_velocity,
        np.zeros((system.shape[0],DATA_SLOTS)) # filler
    ))

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

        return calculate_derivative(
            y.reshape(system.shape), method, theta
        ).ravel()

    solved_state = solve_ivp(
        fun    = differential_equation_functor,
        t_span = (0, end_time),
        y0     = system.ravel(),
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
    feature_positions0 = trajectory[0]

    lines = [
        ax.plot(
            [back[0], front[0]],
            [back[1], front[1]],
            color="black"
        )[0]
        for back, front in zip(
            feature_positions0[:, BACK],
            feature_positions0[:, FRONT]
        )
    ]

    front_scatter = ax.scatter(
        *feature_positions0[:, FRONT].T,
        color="red",
        label="front"
    )

    back_scatter = ax.scatter(
        *feature_positions0[:, BACK].T,
        color="blue",
        label="back"
    )

    ax.set_xlim([-BOUNDS[0], BOUNDS[0]])
    ax.set_ylim([-BOUNDS[1], BOUNDS[1]])
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect("equal")

    def update(frame):
        feature_positions = trajectory[frame]

        fronts = feature_positions[:, FRONT]
        backs = feature_positions[:, BACK]

        ax.set_title(f"t={times[frame]:.2f}")

        for line, back, front in zip(lines, backs, fronts):
            line.set_data(
                [back[0], front[0]],
                [back[1], front[1]]
            )

        front_scatter.set_offsets(fronts)
        back_scatter.set_offsets(backs)

        return [
            front_scatter,
            back_scatter,
            *lines
        ]

    animation = FuncAnimation(
        fig, update, frames=times.size, interval=20, blit=False)

    plt.show()

result = calculate_time_progression(
    system = generate_system(1000, bounds = np.array([50, 50])),
    method = 'brute',
    theta = 0,
    rtol = 1e-2,
    atol = 1e-3,
    integration_method = 'RK23',
    time_step = 0.5,
    end_time = 100,
    show_progress = True
)

feature_trajectory = np.array([
    calculate_feature_positions(system)
    for system in result[0]
])

animate_system(feature_trajectory, result[1])
