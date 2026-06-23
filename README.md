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

### SETUP

All that is required to install FInDS is Docker and GNU Make. For
Windows users, I highly recommend using WSL[1] and following the Linux
instructions to keep your life simple, but this program can be run
purely on windows regardless due to containerization.

DOCKER - LINUX

All that is needed is the Docker CLI. Follow the instructions here:

  <https://docs.docker.com/engine/install/>

DOCKER - WINDOWS AND MACOS

You probably will need Docker Desktop. Follow the instructions here:

  <https://docs.docker.com/desktop/>

GNU MAKE - LINUX AND MACOS

Your system likely comes with Make preinstalled, but if not, follow
the instructions provided by the GNU foundation (ideally, download it
via your package manager):

  <https://www.gnu.org/software/make/>

GNU MAKE - WINDOWS

Alternative packagers like Chocolatey[2] or Scoop[3] for PowerShell
can be used to install Make.


[1] https://learn.microsoft.com/en-us/windows/wsl/install
[2] https://chocolatey.org/
[3] https://scoop.sh/
