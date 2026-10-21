# Complete dotted attachment catalogue through dimension29

The [checker](../scripts/fresh_sqrt5_dotted_rank29_exhaustion_20261015.py)
exhausts the attachment queue of all 32 necessary arithmetic parameter seeds
from the [complete dotted-distance reduction](fresh_sqrt5_dotted_parameters_20261015.md).
Its [certificate](fresh_sqrt5_dotted_rank29_exhaustion_20261015.json) has
**438 marked states and no pending state**. It reduces the dotted part of
the dimension-28/29 argument to **264 explicit face-weight obligations**.
Those weights have not yet been verified; this is not an improved dimension
theorem.

| Processed outcome | States |
| --- | ---: |
| Singular state sealed against further unpaid contact | 250 |
| Singular state with a retained next contact | 4 |
| Nonsingular state with retained next contacts | 156 |
| Sealed nonsingular state excluded by full-star rank | 13 |
| Nonsingular31-root state excluded by ambient rank | 5 |
| Singular31-root state at the dimension boundary | 10 |

The inherited arithmetic upper bound is29. A full quadrilateral star in
dimension d has d+2 roots and Gram rank d+1. Therefore a relevant star has
at most31 roots. At31 roots, a nonsingular Gram matrix already exceeds the
allowed ambient rank; a singular one is a possible zero-contact weight
case in dimension29. This is an explicit restriction to dimensions at most29,
not an arbitrary-rank sealing assertion.

Below that boundary, every possible nonzero next-white contact is generated
in the full finite-angle alphabet, subject only to the reviewed four finite
corners and no-short-path condition. Exact field, physical/conjugate inertia
and cyclic-integrality tests leave422 accepted child occurrences. Every one
has a saved literal marked isomorphism to a queue state. The remaining trials
give3,457 two-negative exclusions,854 finite-cycle field exclusions,334 dotted
cycle field exclusions,66 conjugate exclusions, and13 persistent H4 payments.
The checker also has a nonnegative-root-dependence rejection, but no trial
needed that branch. It makes no negative assertion based on a partial search
for affine subdiagrams.

Deduplication uses exact full-label marked graph isomorphisms within invariant
buckets, retaining the squared dotted distance as part of the identity. It
does not identify geometrically different dotted parameters. There is no
sampling, and the empty queue is checked directly.

In an unpaid actual full star, any nonorthogonal external white component
can be added from its first contact with the core, then along its connected
white edges. Thus it follows a retained branch. A sealed nonsingular core
cannot be completed by orthogonal spectators: their positive-definite Gram
block would leave the whole star nonsingular. A singular state can instead
stop with only orthogonal spectators. Each of the250 sealed singular states,
four nonsealed singular starts, and ten boundary states must therefore receive
the required face payment, with spectator rank limited by the same dimension
bound. Proving those264 obligations pays every previously unresolved dotted
family in the targeted dimensions.

Generation and saved replay pass. Independent review is being performed
separately. All finite-label attachment families and the compatible global
vertex budget remain separate obligations.

```sh
python3 scripts/fresh_sqrt5_dotted_rank29_exhaustion_20261015.py --verify
```
