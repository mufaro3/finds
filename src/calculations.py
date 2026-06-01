from dataclasses import dataclasses
import numpy as np
from numpy.typing import NDArray
from types import SimpleNameSpace
from deepcopy import deepcopy

# define vector3 type
VecF = np.NDArray[np.float64]
    
@dataclass
class Fish:
    """One individual fish model within the school"""
    length: float
    speed: float

    # Doesn't store any derivatives to be compatible with RK4
    dynamic = SimpleNameSpace(
        center_position=None
        head_position=None,
        tail_position=None,
        orientation=None,
        simulation_indx=0
    )

    @staticmethod
    def generate_random_fish_in_box(bounds: VecF) -> Fish:
        # orientation
        theta = np.random.uniform(0, 2 * np.pi)
        cos_phi = np.random.uniform(-1, 1)
        phi = np.arccos(cos_phi)
        
        orientation = np.array([ np.sin(phi) * np.cos(theta),
                                 np.sin(phi) * np.sin(theta),
                                 np.cos(phi) ])

        # position
        position = np.random.random(size=bounds.shape) * bounds

        new_fish = Fish()
        new_fish.dynamic.center_position = position
        new_fish.dynamic.orientation = orientation

        # compute the head and tail positions
        new_fish.dynamic.head_position = new_fish.__calculate_head_pos()
        new_fish.dynamic.tail_position = new_fish.__calculate_tail_pos()

        # we can't compute the velocities just yet
        return new_fish
        
    def __init__(self, position: VecF, orientation: VecF, \
                 length: float = 1, volumetric_flow_rate: float = 4 * np.pi):
        self.length = length
        self.speed = volumetric_flow_rate / (4 * np.pi * length ** 2)

    def __calculate_extrema_pos(self, head: bool) -> VecF:
        """Computes the position of extrema (head/tail)"""
        delta = (self.length * self.orientation) / 2
        if head:
            return self.dynamic.center_position + delta
        else:
            return self.dynamic.center_position - delta
        
    def __calculate_head_pos(self) -> VecF:
        """Calculates the position of the head of the fish."""
        return self.__get_extrema_pos(True)

    def __calculate_tail_pos(self) -> VecF:
        """Calculates the position of the tail of the fish."""
        return self.__get_extrema_pos(False)

    def _calculate_difference_kernel(a: VecF, b: VecF) -> VecF:
        """Computes the difference kernel between two vectors, or 1/r^2"""
        diff = a - b
        return diff / (np.linalg.norm(diff) ** 3)

    def __calculate_interaction_kernel_head(self, other_fish: Fish) -> VecF:
        """Computes the interaction kernel for the head velocity of two fishes"""

        head_pos = self.dynamic.head_position
        other_fish_head_pos = other_fish.dynamic.head_position
        other_fish_tail_pos = other_fish.dynamic.tail_position
        
        head_head = self.__calculate_difference_kernel(head_pos, other_fish_head_pos)
        head_tail = self.__calculate_difference_kernel(head_pos, other_fish_tail_pos)
        
        return head_head - head_tail

    def __calculate_interaction_kernel_tail(self, other_fish: Fish) -> VecF:
        """Computes the interaction kernel for the tail velocity of two fishes"""
        tail_pos = self.dynamic.tail_position
        other_fish_head_pos = other_fish.dynamic.head_position
        other_fish_tail_pos = other_fish.dynamic.tail_position
        
        tail_head = self.__calculate_difference_kernel(tail_pos, other_fish_head_pos)
        tail_tail = self.__calculate_difference_kernel(tail_pos, other_fish_tail_pos)
        
        return tail_head - tail_tail
    
    def __calculate_interaction_kernel_sum(
            self, other_fish_list: list[Fish],
            use_barnes_hut: bool, head: bool) -> VecF:
        """Calculates the interaction kernel sum between a fish
        and a set of its neighbors"""
        interaction_fn = None
        if head:
            interaction_fn = self.__calculate_interaction_kernel_head
        else:
            interaction_fn = self.__calculate_interaction_kernel_tail

        interactions = None
        if use_barnes_hut:
            interactions = None # TODO - Barnes Hut Approximation
        else:
            interactions = \
                np.asarray([ interaction_fn(other_fish) \
                             for other_fish in fish_list ])

        return np.sum(interactions, axis=0)

    def __calculate_extrema_velocity(
            self, other_fish_list: list[Fish],
            use_barnes_hut: bool, head: bool) -> VecF:
        """Computes the velocity of extrema based on the
        orientations and interaction kernels"""
        interaction_kernel_sum = \
            self.__calculate_interaction_kernel_sum(
                other_fish_list, use_barnes_hut, head)
        
        return fish.speed * (fish.orientation + (fish.length ** 2) * \
                             interaction_kernel_sum)
    
    def __calculate_head_velocity(
            self, other_fish_list: list[Fish], use_barnes_hut: bool) -> VecF:
        """Computes the velocity of the head of the fish"""        
        return _calculate_extrema_velocity(other_fish_list, use_barnes_hut, True)

    def __calculate_tail_velocity(
            self, other_fish_list: list[Fish], use_barnes_hut: bool) -> VecF:
        """Computes the velocity of the tail of the fish"""        
        return _calculate_extrema_velocity(other_fish_list, use_barnes_hut, False)

    def __calculate_velocities(self, head_velocity, tail_velocity) -> VecF:
        """Computes the translational and rotational
        velocities of the center-of-mass"""
        
        translational_velocity = (head_velocity + tail_velocity) / 2

        velocity_diff = head_velocity - tail_velocity
        lagrange_multiplier = - 1 / 2 * \
            np.dot(velocity_diff, fish.orientation)
        rotational_velocity = (velocity_diff + 2 * lagrange_multiplier *\
                               self.orientation) / self.length
        
        return translational_velocity, rotational_velocity

    def __compute_derivatives(self, other_fish_list: list[Fish], use_barnes_hut: bool) -> tuple[VecF]:
        """Updates the dynamical state of the fish, with all of the calculations hidden.
        IMPORTANT: This function (and all of its children) must not modify the internal state!"""
        head_velocity = self.__calculate_head_velocity(other_fish_list, use_barnes_hut)
        tail_velocity = self.__calculate_tail_velocity(other_fish_list, use_barnes_hut)

        translational_velocity, rotational_velocity = \
            self.__calculate_velocities(head_velocity, tail_velocity)

        return np.array([translational_velocity, rotational_velocity])

    # Runge-Kutta Simulation Code
    
    def __compute_advanced_state(self, translational_velocity: VecF,
                                 rotational_velocity: VecF, time_step: float) -> Fish:
        """Computes the advanced state of the fish"""
        new_fish_state = deepcopy(self)
        
        new_fish_state.dynamic.center_position += translational_velocity * time_step
        new_fish_state.dynamic.orientation += rotational_velocity * time_step

        new_fish_state.dynamic.head_pos = new_fish_state.__calculate_head_pos()
        new_fish_state.dynamic.tail_pos = new_fish_state.__calculate_tail_pos()

        new_fish_state.dynamic.simulation_index += 1
        
        return new_fish_state
        
    def update_rk4(self, other_fish_list: list[Fish], use_barnes_hut: bool, time_step: float) -> Fish:
        """Computes the next step for the fish according to 4th-Order Runge-Kutta"""
        
        # First Order
        deriv_1 = self.__compute_derivatives(other_fish_list, use_barnes_hut)
        advanced_state_1 = self.__compute_advanced_state(self, *deriv_1, time_step)

        deriv_2 = advanced_state_1.__compute_derivatives(other_fish_list, use_barnes_hut)
        advanced_state_2 = advanced_state_1.__compute_advanced_state(self, *deriv_2, time_step / 2)

        deriv_3 = advanced_state_2.__compute_derivatives(other_fish_list, use_barnes_hut)
        advanced_state_3 = advanced_state_2.__compute_advanced_state(self, *deriv_3, time_step / 2)

        deriv_4 = advanced_state_3.__compute_derivatives(other_fish_list, use_barnes_hut)

        weighted_deriv = (deriv_1 + 2 * deriv_2 + 2 * deriv_3 + deriv_4) / 6
        new_state = self.__compute_advanced_state(self, *weighted_deriv, time_step)

        return new_state
