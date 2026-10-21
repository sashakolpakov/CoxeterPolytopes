# The next singular quadrilateral frontier over Q(√5)

The complete next-white calculation on **51 new singular marked classes**
finds **32 sealed classes paid by the current forest-sensitive rule, 14
sealed classes unpaid by that rule, and five classes with seven surviving
contact profiles**. The 14 unpaid classes supply genuine missing local
weight constraints; they are not excluded simply because the current rule
assigns them zero. No new dimension bound is asserted.

The [checker](../scripts/fresh_sqrt5_singular_extension_frontier_20261015.py)
and [certificate](fresh_sqrt5_singular_extension_frontier_20261015.json)
start with the 61 literal singular first extensions left open by the
[sealed-subcore census](fresh_sqrt5_sealed_subcore_closure_20261015.json).
Every quotient map to the 51 representatives is checked on all labels and
preserves black/white colors and opposite pairs. None of the older frozen
artifacts is changed.

For each representative the physical Gram has inertia (n−2,1,1), and the
conjugate Gram is positive semidefinite with rank n−1. The computation uses
the standalone Fraction Q(√5) arithmetic and finite Dynkin classifier from
the earlier independent contact review. It enumerates every finite-label
profile for a further white node, using the established complete alphabet
2,3,4,5,6,10. Every such contact belongs to an actual spherical corner in a
putative completion, so there is no missing dotted or arbitrary-label
contact in this step.

A short marked path gives persistent baseline payment and is pruned. Every
other partial profile must retain all four positive-definite corner diagrams.
This monotone pruning is valid with the unassigned contacts relaxed to right
angles. The resulting 502 nonzero profiles have the following exact outcomes:

| Outcome | Number |
|---|---:|
| Violates the core relation / conjugate kernel | 377 |
| Contains a certified affine principal diagram | 42 |
| Fails the finite-cycle ground-field condition | 72 |
| Creates an entire H4 payment | 4 |
| Survives all stated necessary tests | 7 |

Affine witnesses are independently checked for positive semidefinite
corank-one Gram and a strictly one-sign kernel. Physical and conjugate
inertias of every final surviving extension are checked exactly. The source
uses a differently scaled spanning tree for each extension; the kernel is
recomputed in that same scaling, so no coefficient is transferred across
incompatible root normalizations.

The kernel rejection also has a direct geometric interpretation. A core
with one negative direction and a Gram kernel has that kernel as a linear
relation among its actual normals: a radical vector in its span is
orthogonal to a timelike direction, so cannot be nonzero null. Every further
normal must respect the relation. Equivalently, its conjugate column must
annihilate the positive-semidefinite core's kernel.

## What sealing proves for arbitrary further whites

For each of the 46 sealed representatives, every nonzero first contact
either fails a necessary condition or produces an already persistent short
path or entire-H4 payment. In a completion not already paid this way, every
external white node must therefore be orthogonal to every core node. A path
through more external whites cannot bypass the argument: its first node
touching the core has one of the enumerated profiles.

Thus the core's marked corner components are entire. The full corner forest
is the displayed core corner forest plus an orthogonal external white
forest. For the present forest-sensitive rule the latter affects payment
only through whether it contains an entire A1 component. The checker tests
both possibilities exactly. This is an all-size sealing implication, not
an extrapolation from one added node or from a finite ambient rank.

The 32 paid sealed classes have total weight at least one in both cases.
The following 14 sealed classes have total weight zero in both cases:

| New index | Connected marked components able to carry payment |
|---|---|
| 15, 18 | two A8 corners |
| 16 | two A7 corners, each accompanied by A1 |
| 28 | two D9 corners |
| 31, 44 | one D9 corner |
| 32, 45 | one B9 corner |
| 33, 46 | two A9 corners |
| 37 | one D8 corner accompanied by A1 |
| 38 | one B8 corner accompanied by A1 |
| 39 | two A8 corners, each accompanied by A1 |
| 43 | one E7+A1 corner and one B8 corner |

The certificate supplies each literal graph, source origin, all four full
corner forests, marked pair orbit, and both weights. Disconnected marked
pairs receive zero under the stated rule. These data can be used to extend
the face-weight constraints. This package does not assert that every
principal star extends to a globally compact arithmetic chamber.

The five nonsealed representatives are indices 2,8,9,29,49. Their seven
profiles remain explicit necessary domains, not successful arbitrary-depth
completions. The 311 nonsingular first extensions and the symbolic dotted
families are outside this calculation.

Generation and saved replay pass, with source manifests, syntax and report
links checked. This is a new primary calculation using previously reviewed
exact helpers; a separate independent replay of this new 51-class stage is
not claimed.

```sh
python3 scripts/fresh_sqrt5_singular_extension_frontier_20261015.py
python3 scripts/fresh_sqrt5_singular_extension_frontier_20261015.py --verify
```
