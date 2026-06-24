Differential Calculation Functions
==================================

This is, essentially, the core module of FINDS. This is where the calculations
that govern the mechanics of the particle simulation occur, and this is written
in accordance with the prior work of Mohamed Niged Mabrouk and Dr. Daniel
Floryan's past work with the far-field model of swimmers
:cite:`mabrouk2024,mabrouk2025`.

The core goal of this module is the computation of the time-derivative of
the system matrix :math:`\mathbf{X}`. Using this calculation, numerical
integration schemes are implemented under the 'simulation' module to allow for
time-stepping.

The following is a simple definition of the far-field model as stated in these
papers alongside descriptions of how this model was implemented in the program
[#]_.


.. [#] Please note that no in-depth justification for this model will be given
       in this documentation! A simple mathematical justification is that the
       leading-order term of the multipole expansion for a body of constant
       volume in fluid is a dipole, but many more justifying reasons are
       stated in the papers of Mabrouk and Floryan
       :cite:`mabrouk2024,mabrouk2025`.

The Far-Field Model of Fish
---------------------------

Each fish-particle consists of two fundamental variables, a center-of-mass
position :math:`\mathbf{x}_c` and an orientation unit vector :math:`\hat{n}`,
and the fish can be represented in full as an array or tuple that stores both
of these variables, :math:`\mathbf{s} = (\mathbf{x}_c, \hat{n})`, and the
system of swimmers :math:`\mathbf{X}` can be defined as a matrix consisting of
each of these fish vectors as rows. Thus, for a system of size :math:`N`,
the system matrix is :math:`N \times 6` or of shape :math:`(N,6)`.

A simple analogy for the far-field model for fish is as treating swimmers like
bar-magnets in free-space, or simple linear dipoles separated by the fish
length :math:`\ell` such that the head (or front) of the fish is a source (or
the north end of the magnet) given that the head pushes fluid away radially
and the tail (or back) of the fish is a sink (or south end of the magnet) given
that water re-accumulates at the tail of the fish.

Mathematically, these two ends can differ only up to a sign (so the head and
tail also form positive and negative ends respectively) due to the conservation
of mass of the fluid (the two ends must cancel each other out). Due to this
duality, the sources and sinks can be considered more generally as "features"
of a fish.

Computing the Derivative
------------------------

For a given fish :math:`\mathbf{s}`, the features (front/source and back/sink)
can be computed as

.. math::
   \begin{aligned}
     \mathbf{x}_f &= \mathbf{x}_c + \vec{\delta} \\
     \mathbf{x}_b &= \mathbf{x}_c - \vec{\delta}
   \end{aligned}

where

.. math::
     \vec{\delta} = \frac{1}{2} \ell \hat{n}

The core equation governing the motion of bodies within the simulation is the
induced velocity for a particle a displacement :math:`\mathbf{r}` away from
a feature (some source or sink of a fish),

.. math::
   \mathbf{u} = \pm \frac{\sigma}{4\pi} \frac{\mathbf{r}}{|\mathbf{r}|^3}

where :math:`\sigma` is the volumetric flow rate of the source or sink.
Sources produce positive velocities as they repel bodies outward and sinks
produce negative velocities as they attract bodies toward them.

In essence, this equation is used to develop the differential time-derivative
of the full state, as the velocity for a single feature is defined as the
total velocity contribution from all other features in the system, and the
feature velocity and orientation of a fish define the velocity of its
orientation.

This idea leads to the feature velocity equation

.. math::
   \mathbf{v}_\alpha = U \hat{n} + \frac{\sigma}{4\pi} \mathbf{h}_\alpha

that computes the velocity of some feature :math:`\alpha` (which would be
either the front or back) for a fish. This is comprised of two terms, an
internal contribution to the velocity resulting from the interaction between
this feature and its polar opposite on the same fish (also known as the
self-propelled velocity) defined as

.. math::
   U \hat{n} = \frac{\sigma}{4 \pi \ell^2} \hat{n},

and the external contribution term that encodes the total velocity
contribution from all features (sources and sinks) external to this fish,
where the vectors :math:`\mathbf{h}` and :math:`\mathbf{c}` are defined to
represent the value of the original term
:math:`\frac{\mathbf{r}}{|\mathbf{r}|^3}`, which I call the "interaction"
between features.

In this context, :math:`\mathbf{h}_{\alpha,i}` refers to the total interaction
between feature :math:`\alpha` of fish :math:`i` with all of the other features
(barring its polar opposite on fish :math:`i`) within the greater system
:math:`\mathbf{X}`,

.. math::
   \frac{\sigma}{4\pi} \mathbf{h}_{\alpha,i} =
   \sum_{j \neq i}^N \frac{\sigma}{4\pi} \mathbf{h}_{\alpha,ij}.

Furthermore, :math:`\mathbf{h}_{\alpha,ij}` refers to the total interaction
between feature :math:`\alpha` of fish :math:`i` and both features of fish
:math:`j` (both the front and the back), and this is always defined as the
individual interaction to the front of fish :math:`j` minus the individual
interaction to the back of fish :math:`j` due to the aforementioned sign
change between sources and sinks,

.. math::
   \mathbf{h}_{\alpha,ij} = \mathbf{c}_{\alpha f} - \mathbf{c}_{\alpha b}.

Lastly, :math:`\mathbf{c}_{\alpha \beta}` is defined as the individual
interaction between feature :math:`\alpha` of fish :math:`i` and feature
:math:`\beta` of fish :math:`j` within :math:`\mathbf{X}`, and this is
defined simply as

.. math::
   \mathbf{c}_{\alpha\beta} =
   \frac{\mathbf{r}_{\alpha\beta}}{|\mathbf{r}_{\alpha\beta}|^3}

where :math:`\mathbf{r}` is defined as the displacement between the two
features

.. math::
   \mathbf{r}_{\alpha\beta} = \mathbf{x}_{\alpha,i} - \mathbf{x}_{\beta,j}.

With this, the derivative for one fish can be defined as the derivative of
the center-of-mass position and the orientation

.. math::
   \dot{\mathbf{s}} = \left(\dot{\mathbf{x}}_c, \frac{d \hat{n}}{dt}\right),

and equivalently, :math:`\dot{\mathbf{x}}_c` is the translational velocity and
:math:`\frac{d \hat{n}}{dt}` is the rotational velocity of the fish. The
translational velocity for a fish is just the average of the front and back
velocities,

.. math::
   \dot{\mathbf{x}}_c = \frac{\mathbf{v}_f + \mathbf{v}_b}{2}

and the rotational velocity is defined as

.. math::
   \frac{d \hat{n}}{dt} = \frac{\delta \mathbf{v} + 2 \lambda \hat{n}}{\ell}

where :math:`\Delta v = \mathbf{v}_f - \mathbf{v}_b` is the difference in
velocity between the front and back of the fish and

.. math::
   \lambda = \frac{-\Delta \mathbf{v} \cdot \hat{n}}{2}

is a Lagrange multiplier used to ensure that the length of the fish
:math:`\ell` is kept constant.

Functions
---------

.. automodule:: finds.calculations
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: __init__
