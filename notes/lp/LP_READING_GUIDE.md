# Reading the linear programming proof

The [36-page manuscript](../../arxiv/linear_programming_dimension_bounds/main.pdf)
uses linear programming to combine previously proved local inequalities.
It treats finite-volume, compact and compact arithmetic Coxeter polytopes.
The recovered dimension bounds are **57, 17 and 15**, respectively.
The [submission guide](../../arxiv/linear_programming_dimension_bounds/README.md)
identifies the current source, PDF and package.

## Start with one compact configuration

Fix a vertex diagram and the additional facet conditions required by the
local theorems. Choose a vertex farthest from a fixed interior point.
The geometric construction gives a tangent vector $x$ with
$\|x\|^2<1$. The connected components of the vertex diagram span
mutually orthogonal subspaces $W_j$. Put

$$
e_j=\|\operatorname{proj}_{W_j}x\|^2.
$$

These nonnegative scalars are the LP variables. Every inequality uses the
same vector $x$, and $\sum_j e_j=\|x\|^2<1$.
The local theorems give inequalities $Ae\ge c$: a row records which
component contributions are summed and their rational lower bound.
Consider

$$
m=\min\{\mathbf1^Te:e\ge0,\ Ae\ge c\}.
$$

If $m\ge1$, the configuration cannot occur. If $m<1$, these
inequalities alone do not exclude it. A feasible LP vector need not have a
geometric realization.

## How the dual proves the contradiction

The dual seeks nonnegative multipliers $\alpha$ satisfying

$$
A^T\alpha\le\mathbf1,
\qquad c^T\alpha\ge1.
$$

The first condition limits the coefficient of each $e_j$ to one in
the chosen combination. It prevents counting an overlapping component
twice. Multiplication gives

$$
\mathbf1^Te\ge\alpha^TAe\ge\alpha^Tc\ge1,
$$

contradicting the geometric strict upper bound. Exact rational checks
verify the multipliers; a numerical solver's status is insufficient.

The manuscript's homogeneous formulation combines the configuration
classes into one feasibility problem. Its scaling variables express the
strict inequality without imposing a fixed positive tolerance. Local
geometric exclusions set the corresponding scaling coefficient to zero.

## Finite volume and the meaning of a mark

Finite volume uses incidence counts rather than projection norms. A
weight function is summed over the corners of each two-dimensional face
and over the corners at each vertex. The LP seeks lower bounds on the
face sums and upper bounds on the vertex sums that contradict the
incidence inequality.

A corner, or planar angle, is a pair $(v,F)$. At a finite vertex,
the diagram nodes represent the incident facets. The two **marked nodes**
represent the facets not containing $F$; the other $n-2$ facets
contain it. Marks identify the chosen incidence. An **edge weight** instead
comes from the Gram matrix and describes two facet hyperplanes. Marked
nodes need not be adjacent. Marking a component in the compact case
similarly records an exceptional configuration, rather than an edge weight.

## From configurations to dimension bounds

An exhaustive family must contain every configuration arising from a
hypothetical polytope. Arithmetic reductions and the local geometric
inequalities remain proved inputs to the current LP argument. The LP
combines them; it does not rederive them from elementary Gram conditions.

Monotonicity requires one fixed family labelled by ambient dimension.
Raising the threshold removes classes, so a feasible point extends by zero
to the lower threshold. Thus infeasibility excludes every higher threshold
within the covered range. Finite volume includes analytic estimates in all
higher dimensions. The compact programs cover dimensions through 29 and
use Vinberg's theorem above that range. See the
[three-case formulation](UNIFORM_LP_PARADIGM.md) and
[exact recovery records](existing_proof_lp_recovery_results_20261118.md).
