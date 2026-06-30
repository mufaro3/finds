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

The fundamental approximation of the far-field model is treating fish as
self-propelled dipoles moving in free space, with the head of the fish being
the source and the tail of the fish being the sink. The head and the tail for
each fish are separated :math:`\ell` apart, or the length of the fish.

Mathematically, these two ends can differ only up to a sign (so the head and
tail also form positive and negative ends respectively) due to the conservation
of mass of the fluid (the two ends must cancel each other out). Due to this
duality, the sources and sinks can be considered more generally as "features"
of a fish.

We store all of the fish in the system in a matrix :math:`\mathbf{X}` where
each fish vector comprises a row, and this is called the system (or state)
matrix.

Computing the Derivative of the System Matrix
---------------------------------------------

For a given fish :math:`\mathbf{X}_i`, the features (front/source and back/sink)
can be computed as

.. math::
   :label: feature_positions

   \mathbf{x}_{f,i} &= \mathbf{x}_{c,i} + \vec{\delta}_i \\
   \mathbf{x}_b &= \mathbf{x}_{c,i} - \vec{\delta}_i

where

.. math::
   :label: fish_delta

   \vec{\delta}_i = \frac{1}{2} \ell \hat{n}_i

and this is returned as the following matrix

.. math::
   :label: feature_positions_matrix

   \mathbf{F} = \begin{bmatrix}
   x_{f1} & y_{f1} & z_{f1} & x_{b1} & y_{b1} & z_{b1} \\
   x_{f2} & y_{f2} & z_{f2} & x_{b2} & y_{b2} & z_{b2} \\
   \vdots & \vdots & \vdots & \vdots & \vdots & \vdots \\
   x_{fN} & y_{fN} & z_{fN} & x_{bN} & y_{bN} & z_{bN}
   \end{bmatrix}

representing the feature positions.

The core equation governing the motion of bodies within the simulation is the
induced velocity for a particle a displacement :math:`\mathbf{r}` away from
a feature (some source or sink of a fish),

.. math::
   :label: induced_velocity

   \mathbf{u} = \pm \frac{\sigma}{4\pi} \frac{\mathbf{r}}{|\mathbf{r}|^3}

where :math:`\sigma` is the volumetric flow rate of the source or sink.
Sources produce outward velocities as they repel bodies outward and sinks
produce inward velocities as they attract bodies toward them.

In essence, this equation is used to develop the differential time-derivative
of the full state, as the velocity for a single feature is defined as the
total velocity contribution from all other features in the system, and the
feature velocity and orientation of a fish define the velocity of its
orientation.

This idea leads to the feature velocity equation

.. math::
   :label: feature_velocity

   \mathbf{v}_{\alpha,i} = U \hat{n}_i + \frac{\sigma}{4\pi}
   \mathbf{h}_{\alpha,i}

(see :py:func:`finds.calculations.calculate_feature_velocities`)
that computes the velocity of some feature :math:`\alpha` (which would be
either the front or back) for a fish. This is comprised of two terms, an
internal contribution to the velocity resulting from the interaction between
this feature and its polar opposite on the same fish (also known as the
self-propelled velocity) defined as

.. math::
   :label: interal_contribution

   U \hat{n}_i = \frac{\sigma}{4 \pi \ell^2} \hat{n}_i,

and the external contribution term that encodes the total velocity
contribution from all features (sources and sinks) external to this fish,
where the vectors :math:`\mathbf{h}` and :math:`\mathbf{c}` are defined to
represent the value of the original term :math:`\mathbf{r}/r^3`, which I call
the "interaction" between features.

In this context, :math:`\mathbf{h}_{\alpha,i}` refers to the total interaction
between feature :math:`\alpha` of fish :math:`i` with all of the other features
(barring its polar opposite on fish :math:`i`) within the greater system
:math:`\mathbf{X}`,

.. math::
   :label: external_contribution

   \frac{\sigma}{4\pi} \mathbf{h}_{\alpha,i} =
   \sum_{j \neq i}^N \frac{\sigma}{4\pi} \mathbf{h}_{\alpha,ij}.

The sum :math:`\sum_{j \neq i}^N \dots` is the core bottleneck of this program,
as this computation for each fish produces an :math:`\mathcal{O}(N^2)` runtime
when fully pairwise (see
:py:func:`finds.calculations.compute_interaction_pairwise`). Thus, this
is where the system diverges into using a tree-based computation method, the
Barnes-Hut algorithm (see
:py:func:`finds.calculations.compute_interaction_barnes_hut` and
:ref:`barnes_hut_section`).

Furthermore, :math:`\mathbf{h}_{\alpha,ij}` refers to the total interaction
between feature :math:`\alpha` of fish :math:`i` and both features of fish
:math:`j` (both the front and the back), and this is always defined as the
individual interaction to the front of fish :math:`j` minus the individual
interaction to the back of fish :math:`j` due to the aforementioned sign
change between sources and sinks,

.. math::
   :label: feature_interaction

   \mathbf{h}_{\alpha,ij} = \mathbf{c}_{\alpha f} - \mathbf{c}_{\alpha b}.

(see :py:func:`finds.calculations.calculate_fish_interaction`)
Lastly, :math:`\mathbf{c}_{\alpha \beta}` is defined as the individual
interaction between feature :math:`\alpha` of fish :math:`i` and feature
:math:`\beta` of fish :math:`j` within :math:`\mathbf{X}`, and this is
defined simply as

.. math::
   :label: interaction_vector

   \mathbf{c}_{\alpha\beta} =
   \frac{\mathbf{r}_{\alpha\beta}}{|\mathbf{r}_{\alpha\beta}|^3}

(see :py:func:`finds.calculations.calculate_feature_interaction`)
where :math:`\mathbf{r}` is defined as the displacement between the two
features

.. math::
   :label: feature_displacement_vector

   \mathbf{r}_{\alpha\beta} = \mathbf{x}_{\alpha,i} - \mathbf{x}_{\beta,j}.

With this, the derivative for one fish can be defined as the derivative of
the center-of-mass position and the orientation

.. math::
   :label: fish_derivative

   \dot{\mathbf{X}}_i = \left(\dot{\mathbf{x}}_{ci}, \frac{d \hat{n_i}}{dt}\right),

(see :py:func:`finds.calculations.calculate_system_derivative`)
and equivalently, :math:`\dot{\mathbf{x}}_c` is the translational velocity and
:math:`\frac{d \hat{n}}{dt}` is the rotational velocity of the fish. The
translational velocity for a fish is just the average of the front and back
velocities,

.. math::
   :label: fish_translational_derivative

   \dot{\mathbf{x}}_c = \frac{\mathbf{v}_f + \mathbf{v}_b}{2}

and the rotational velocity is defined as

.. math::
   :label: fish_rotational_derivative

   \frac{d \hat{n}}{dt} = \frac{\Delta \mathbf{v} + 2 \lambda \hat{n}}{\ell}

where :math:`\Delta \mathbf{v} = \mathbf{v}_f - \mathbf{v}_b` is the difference in
velocity between the front and back of the fish and

.. math::
   :label: lagrange_multiplier

   \lambda = \frac{-\Delta \mathbf{v} \cdot \hat{n}}{2}

is a Lagrange multiplier used to ensure that the length of the fish
:math:`\ell` is kept constant.

.. _barnes_hut_section:

Barnes-Hut Clustering Approximation
-----------------------------------

The Barnes-Hut algorithm simplifies calculating the interaction vectors by
clustering fish that are sufficiently "far away" from a given fish (as
determined by the Barnes-Hut ratio :math:`\theta`). A complete, interactive
description of the Barnes-Hut algorithm can be found online
:cite:`heer_barnes_hut`.

This begins by constructing an Octree that partitions the three-dimensional
space around the origin. The fish-particles are inserted in list order, and for
each additional point, the Octree expands by further subdividing the three-
dimensional space. Then, once the Octree is fully built, each node of the tree
(essentially representing every possible division of the three-dimensional
space) is "clustered," meaning that several of the fish-particles are computed
into a fish-particle representing the average of all of them.

For example, if we have fish :math:`\mathbf{X}_i` and fish :math:`\mathbf{X}_j`
in :math:`\mathbf{X}`, their clustered form would simply be the average of the
two, :math:`(\mathbf{X}_i + \mathbf{X}_j)/2`. Then, we traverse this tree for
each fish. If the node is a leaf, then we automatically compute the
interaction. For each branch node in the octree, we then compute a size-to-
distance ratio :math:`\phi`. This is computed as the fraction of the "size" or
the side length of the cube comprising the subdivision over the distance from
the center-of-mass of the subdivision (or the position of the cluster).

For example, if we have a subdivision storing fish 1 through :math:`k`
clustered into a cluster-fish with position :math:`\mathbf{x}_{c}` with
a side length of size :math:`l`, then the Barnes-Hut ratio with respect
to a fish at position :math:`\mathbf{x}_{c0}` would be

.. math::
   :label: barnes_hut_ratio

   \phi = \frac{l}{||\mathbf{x}_c - \mathbf{x}_{c0}||}.

Once we've calculated :math:`\phi` for the given node, if :math:`\phi \ge
\theta`, then we continue travering the tree to the children of the node. If
:math:`\phi < \theta`, then we compute the interaction.

Functions
---------

.. automodule:: finds.calculations
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: __init__
