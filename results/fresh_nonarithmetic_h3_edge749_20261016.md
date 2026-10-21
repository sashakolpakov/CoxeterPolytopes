# Recentered endpoint and coupled-facet search near edge [7,49]

This fresh numerical experiment rebuilds the actual H3 geometry near a
small exterior region. It finds compatible two-cap proposals whose first
new facet is a finite triangle, but none cuts more than one original open
edge end. The whole H3 remains open. No proposal is promoted to an exact
cap, compact H3 or finite-volume ambient chamber.

## Why the chart changed

The preceding [seven-profile search](fresh_nonarithmetic_h3_a2_profiles_20261016.md)
left one end of old edge `[7,49]` uncovered after finite-incidence pruning.
It also discarded 3,023 sampled occurrences with uncertain floating-point
norms. In the old Lorentz chart, the relevant normals could have coordinates
of size `10^4` or larger. That finite-bank deficit was not an obstruction to
the continuous geometry.

The [recentered source](../scripts/fresh_nonarithmetic_h3_edge749_20261016.py)
and [saved data](fresh_nonarithmetic_h3_edge749_20261016.json) reconstruct all
55 old physical normals from the frozen exact sources and recompute the
20 face rays with high-precision linear algebra. A timelike chart vector
is chosen near the midpoint of the two ideal endpoints of edge `[7,49]`,
with a small displacement toward the previous interior chart vector.
Its products with every non-S6 old wall are checked negative.

The common center is also changed. Write the old center as its white
S6 component `H_white` plus its orthogonal component. The new center is

`X = H_white + k*T_new`,

where `k` is chosen large enough that `X` is timelike and strictly inside
every old wall. These signs are checked at high precision. For a new A2
projection `p`, its residual time-coordinate test is
`y_0 > (p,X)/k`; every native old pairing still includes `(p,n_j)`.
These are numerical common-center checks, not rational interval
certificates. Retained walls are unchanged, but retaining the original
particular center is no longer imposed.

The same bounded endpoint calculation is applied to all seven profiles,
including the zero profile. It checks 7,371 recipes, 732 real endpoint
branches and 65,880 isolated third-contact equations. There are 402
admissible open families and 4,823 numerical candidates. The original
[exact profile lemma](fresh_nonarithmetic_h3_a2_profiles_20261016.md)
is unchanged. The chart, center and sample parameters have changed, so this
is not a one-to-one precision replay of the preceding candidate bank.

## Finite-incidence screen

The [joint source](../scripts/fresh_nonarithmetic_h3_edge749_20261016_joint.py)
and [joint data](fresh_nonarithmetic_h3_edge749_20261016_joint.json) reapply
native Coxeter compatibility and projected finite-neighbor tests. They
skip 268 norm-uncertain occurrences and deduplicate one row, leaving 4,554
candidates. The finite new-new graph has 28 edges. Iterated polygon-degree
pruning leaves 1,776 candidates.

The earlier target on `[7,49]` has surviving cutters. However, none of the
remaining candidates cuts ideal-end target0 on old edge `[7,30]`. The
373,930-row finite compatible-cover model is infeasible. This is still a
bounded floating-point result, not an all-parameter or fixed-wall theorem.

## A coordinated continuation at the missing end

The old support `S6+[7,30]` is `H4+H4`. A finite endpoint there must therefore
be supplied by a wall right to those eight roots: the positive inverse of
either H4 makes any nonright contact exceed the available norm cost.
Seven sampled endpoint caps `r` are present in the recentered bank, but
they have only two old finite intrinsic neighbors.

The [closure source](../scripts/fresh_nonarithmetic_h3_edge749_20261016_closure.py)
and [closure data](fresh_nonarithmetic_h3_edge749_20261016_closure.json) fix
each such `r` and solve a second cap `s` with prescribed finite contacts to
`r,7,30`. All seven A2 profiles and labels 2 through 10 are included.
The two additional corners of the proposed triangular `r`-facet must have
positive-definite residual Grams. The norm equation then has at most two
branches, both tested against the common center and all 55 old native gates.

Of 35,721 recipes, 35,224 fail a necessary finite-corner test; 63 have no
real branch. The remaining 868 branches produce 323 compatible numerical
pairs. In 261 of them, clipping gives exactly three finite vertices on
`r` and no exterior ray on that facet. The other 62 have near-null output
and remain unresolved. The 261 finite cases give a whole-H3 census of
16 finite vertices and **eight exterior rays**, compared with seven
exterior rays initially.

Every accepted pair cuts only the same initial edge-end target0. Thus this
is a genuine coupled finite-facet mechanism, but it has not closed a larger
exterior region. Exactifying an arbitrary one of these local triangles
would not establish the construction objective. The next coordinated
target is the open edges on the replacement `s`-facet, allowing all seven
A2 profiles and retaining the native offsets. A polygon having three
possible finite neighbors is not itself a compactness certificate.

In particular, the H4+H4 argument forcing `r` right to its endpoint support
does **not** prohibit a third side `s` of the triangular `r`-facet. Its
finite corners require `PD(S6+r+7+s)` and `PD(S6+r+30+s)`, not
`PD(S6+7+30+s)`. The source tests the correct two corner conditions.

## Replay and scope

Generation and fresh saved replay pass for all three sources. Syntax,
source manifests and report links are checked separately. They reuse the
frozen profile and numerical clipping helpers; no frozen input is edited.
All claims about candidate acceptance, finite corners and binary-program
failure in this report are numerical. An exact cap certificate, complete
continuous-family coverage and a finite-volume example remain open.

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_nonarithmetic_h3_edge749_20261016.py --verify
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_nonarithmetic_h3_edge749_20261016_joint.py --verify
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_nonarithmetic_h3_edge749_20261016_closure.py --verify
```
