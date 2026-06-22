# FINDS: The Fish INteraction and Dynamics Simulator

FINDS was written to solve the problem of the unscalability of fish-like
particle simulations used to understand the mechanics of fish swimming in
schools, birds flying in groups, drones/bugs flying in swarms, and so on.

This fundamentally follows the N-body problem, as each fish-like particle
(which I have equivalently called "fish automatons," "autofish," "fish-
particles," and even just "particles" in various places) affects each other
fish-like particle in the simulation. At a group size of N=3, an analytical
solution is no longer possible, and numerical calculation is needed. As the
group size N increases, the computational complexity in performing these
numerical calculations grows as a function of N^2.

As such, the fundamental goal of 'Fish' was to simplify the calculations
to reduce the effect of this bottleneck on Floryan's research through the
implementation of approximation algorithms and speed tricks commonly used
in Particle dynamics, starting with the Barnes-Hut approximation, which
decreases the computational complexity to O(N log N).

By August 20 onward, future updates will be unlikely to continue to this
repository.
