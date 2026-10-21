# Distinct arithmetic cap profiles cannot always coexist

The [exact certificate](fresh_sqrt5_cap_profile_cliques_20261017.json)
and [producer](../scripts/fresh_sqrt5_cap_profile_cliques_20261017.py)
add36 capacity-one clique inequalities to the preceding815 cap columns,
giving851 columns in total. Four new inequalities violate the original
full-model dual. These constraints concern different cap profiles at
one actual simplex type, beyond repetitions or symmetries of one profile.

## Geometric and arithmetic argument

Start with the capacity-one profiles from the
[cross-orientation package](fresh_sqrt5_cross_orientation_cap_capacity_20261017.md).
Group them by the exact base type, retaining its canonical core and
finite spectators. Compare every orientation of one profile with every
orientation of the other, under a single fixed identification of this
core. For each selected pair, their base-vertex subsets intersect.
Thus they cannot occur at the same isolated A1 coface: the inherited
vertex/coface flag capacity would be exceeded.

At distinct isolated A1 cofaces, their cap-wall product is the product
of their core projections, since the two coface vectors are orthogonal
to each other and to the core. Its square fails arithmetic two-cycle
integrality in every selected orientation pair. This excludes coexistence
at distinct cofaces also. The dotted cap/coface parameters of the two
profiles may differ; those terms disappear from this cross product.

Each individual profile already has capacity one. A collection of
pairwise incompatible profiles therefore has total capacity one. Adding
their incidence coefficients gives a valid signed column with one unit
of the common containing-simplex resource. The producer checks disjoint
source-base lists before adding coefficients. No arbitrary scaling or
independence between separate clique inequalities is assumed.

## Exact selected domain and output

There are80 pairs of capacity-one profiles sharing a base type. We
select cores whose physical normalized Gram matrix lies in Q(sqrt(5))
with labels3 and5, and cap contacts with labels2 and3. Eight pairs with
dotted cores lie outside this selected calculation. Two other pairs have
an integral cross-product square and are not excluded by this test.
The remaining70 pairs are incompatible by the argument above. The data
record every tested pair and every skipped or non-excluded case.

The checker reconstructs the exact core inverse, multiplies it against
the Gram matrix, verifies both individual projection norms, and records
the cross-product square, trace and norm. Integral-basis and trace/norm
tests agree. All70 excluded pairs pass the arithmetic obstruction;
no finite-angle alphabet argument is used in this package.

The resulting graph has36 maximal cliques of size at least two:22 pairs,
12 triples and2 quadruples. Maximality is checked exhaustively within
each finite saved profile group; this is not a classification of all
possible geometric cap profiles. Thirty-four new columns have positive
load at the original dual and four violate it. Together with the42
violated predecessor columns,46 columns are violated in this851-column
family. This count concerns the unchanged diagnostic masses.

Generation and fresh saved replay pass, with no optimizer invoked.
The producer shares the old reviewed exact field helper and binds that
source explicitly. An independent family review and any repair at changed
face masses have separate scopes. This result supplies constraints, not
a new optimum, a finite-volume example or an ambient dimension bound.
