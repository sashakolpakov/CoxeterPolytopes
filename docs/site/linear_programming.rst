Linear programming and the dimension bounds
===========================================

The :download:`36-page LP manuscript <../../arxiv/linear_programming_dimension_bounds/main.pdf>`
combines proved local inequalities to recover three existing bounds.

.. list-table:: Bounds recovered by the current LP arguments
   :header-rows: 1
   :widths: 55 20 25

   * - Coxeter polytopes
     - Dimension bound
     - Excluded dimensions
   * - Finite volume
     - 57
     - 58 and higher
   * - Compact
     - 17
     - 18 and higher
   * - Compact arithmetic
     - 15
     - 16 and higher

The LP supplies the global deduction. Local geometry, arithmetic reductions
and completeness of the configuration lists remain inputs.

One compact configuration
-------------------------

Fix a vertex diagram and the additional facet conditions required by the
local theorems. Choose a vertex farthest from a fixed interior point. The
construction gives a tangent vector :math:`x` with :math:`\|x\|^2<1`.
The connected components of the vertex diagram span mutually orthogonal
subspaces :math:`W_j`. Set

.. math::

   e_j=\|\operatorname{proj}_{W_j}x\|^2.

These nonnegative scalars are the LP variables. Every inequality uses
the same original vector :math:`x`, and
:math:`\sum_j e_j=\|x\|^2<1`.
The local inequalities have the form :math:`Ae\ge c`. Each row specifies
which component contributions are summed and their rational lower bound.
The primal problem is

.. math::

   m=\min\{\mathbf 1^Te:e\ge0,\ Ae\ge c\}.

Infeasibility or :math:`m\ge1` excludes this configuration. A feasible
vector with total below one does not establish geometric realizability.

The dual certificate
--------------------

Find nonnegative multipliers :math:`\alpha` such that

.. math::

   A^T\alpha\le\mathbf1,
   \qquad c^T\alpha\ge1.

The first condition keeps each component's coefficient at most one,
including when several inequalities involve it. Consequently,

.. math::

   \mathbf1^Te\ge\alpha^TAe\ge\alpha^Tc\ge1,

contradicting the geometric upper bound. Exact rational checks verify
these multipliers. The manuscript combines configuration classes into a
homogeneous feasibility problem; scaling expresses the strict upper bound
without choosing a fixed numerical tolerance.

Finite volume and diagram marks
--------------------------------

The finite-volume LP instead assigns weights to vertex–two-face
incidences. The same weight function is summed over each two-dimensional
face and over the incidences at each vertex. Lower bounds for the face
sums and upper bounds for the vertex sums contradict an incidence
inequality.

An incidence, or planar angle, is a pair :math:`(v,F)`. At a finite vertex,
the diagram nodes represent the incident facets. The two marked nodes
represent the facets not containing :math:`F`; the other :math:`n-2`
facets contain it. A mark selects nodes. An edge weight describes the
relative position of two facet hyperplanes through their Gram entry.
Marked nodes need not be adjacent.

Coverage and monotonicity
-------------------------

Every hypothetical polytope must produce a configuration in the stated
domain and satisfy its inequalities. The current LP does not derive
these local theorems from elementary Gram conditions.

For monotonicity, fix one configuration family labelled by ambient
dimension. Raising the dimension threshold removes classes. A feasible
point at the higher threshold extends by zero to the lower threshold;
infeasibility therefore excludes all higher thresholds in the covered
range. Finite volume includes analytic estimates for the higher-dimensional
tail. The compact programs cover dimensions through 29 and use Vinberg's
theorem above that range.

See the :download:`short reading guide <../../notes/lp/LP_READING_GUIDE.md>`,
:download:`three-case formulation <../../notes/lp/UNIFORM_LP_PARADIGM.md>` and
:download:`exact recovery records <../../notes/lp/existing_proof_lp_recovery_results_20261118.md>`
for the definitions, retained inputs and execution evidence.
