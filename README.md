# FINDS: The Fish INteraction and Dynamics Simulator

FINDS was written to solve the problem of the unscalability of fish-like
particle simulations used to understand the mechanics of fish swimming in
schools, birds flying in groups, drones/bugs flying in swarms, and so on.

This fundamentally follows the N-body problem, as each fish-like particle
(which I have equivalently called "fish automatons," "autofish," "fish-
particles," and even just "particles" in various places) affects each other
fish-like particle in the simulation. At a group size of $N=3$, an analytical
solution is no longer possible, and numerical calculation is needed. As the
group size $N$ increases, the computational complexity in performing these
numerical calculations grows as a function of $N^2$.

As such, the fundamental goal of 'Fish' was to simplify the calculations
to reduce the effect of this bottleneck on Floryan's research through the
implementation of approximation algorithms and speed tricks commonly used
in Particle dynamics, starting with the Barnes-Hut approximation, which
decreases the computational complexity to $O(N \log N)$.

Long term, the goal for the project is to try and reach computational parity
with modern particle dynamics codes used in Astrophysics (and the like), such
as [REBOUND](https://rebound.hanno-rein.de/),
[GADGET-4](https://wwwmpa.mpa-garching.mpg.de/gadget4/), and
[similar](https://github.com/pmocz/awesome-astrophysical-simulation-codes)
applied to fish simulation. This is, of course, an ambitious goal for a
one-man project, but we'll see how this fares.

## SETUP

All that is required to install FInDS is Docker and GNU Make. For
Windows users, I highly recommend using
[WSL](https://learn.microsoft.com/en-us/windows/wsl/install) and following the
Linux instructions to keep your life simple, but this program can be run
purely on windows regardless due to containerization.

### Setup Example for Debian-Derived Systems

1. Download prerequisites (docker and make)

```shell
sudo apt update
sudo apt install git make
sudo apt remove $(dpkg --get-selections docker.io docker-compose docker-doc \
	podman-docker containerd runc | cut -f1)

# Add Docker's official GPG key:
sudo apt install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/debian/gpg -o \
	/etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
sudo tee /etc/apt/sources.list.d/docker.sources <<EOF
Types: deb
URIs: https://download.docker.com/linux/debian
Suites: $(. /etc/os-release && echo "$VERSION_CODENAME")
Components: stable
Architectures: $(dpkg --print-architecture)
Signed-By: /etc/apt/keyrings/docker.asc
EOF

sudo apt update
sudo apt install docker-ce docker-ce-cli containerd.io docker-buildx-plugin \
	docker-compose-plugin
```

2. Clone the Repository

```shell
git clone https://github.com/mufaro3/finds.git
```

3. `cd` into the repository and build the container and the documentation

```shell
cd finds
make build
make docs
```

4. FINDS is now ready for use. Take a look at `make help` for options on what
   you can do further from there.

### Docker

#### Docker - Linux

All that is needed is the Docker CLI. Follow the instructions on the
[Docker Website](https://docs.docker.com/engine/install/).

#### Docker - Windows and MacOS

You probably will need Docker Desktop. Follow the instructions on the
[Docker Website](https://docs.docker.com/desktop/).

### GNU Make

#### GNU Make - Linux and MacOS

Your system likely comes with Make preinstalled, but if not, follow
the instructions provided by the GNU foundation (ideally, download it
via your package manager) at [GNU](https://www.gnu.org/software/make/).

## Acknowledgements

Thanks to Mohamed Niged Mabrouk, Professor Daniel Floryan, and Alvin Ng for
their guidance on this project. In particular, Alvin Ng's
[grav_sim](https://github.com/alvinng4/grav_sim) was used extensively as a
reference work for developing FINDS due to its similarity in implementation.
