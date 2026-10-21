# 393 new silver walls at `x=4+3sqrt(2)`

The intrinsic norm-two Vinberg search on the designed silver K3 carrier is
now exhausted through

\[
x=4+3\sqrt2.
\]

This shell adds exactly 393 walls relative to the certified 333-wall prefix,
giving a 726-wall acute prefix.  New wall 622 is a second finite endpoint of
the targeted `A2+F4+B2` edge, whose finite endpoint set is now `(33,622)`.
This is a local edge closure only: the 726-wall list is not known to be the
complete simple-root system and no global Vinberg stopping or volume theorem
follows.

The frozen manifest is
`data/candidates/silver_k3_norm2_vinberg_extension_4_plus_3sqrt2.json`; its
canonical digest is

```text
391b00b000fb57e8073c873a8db0685d42138dd5cf795f3150bfe6ab430b62c1
```

The block-factorized replay has digest

```text
87313bb8f9911ba9a6102c622475cfe09a6940685cafe32f42ba57f0ab449adf
```

## Exact shell enumeration

Write `s=sqrt(2)`.  The finite-weight equation is

\[
z^T K^{-1}z=2+(1+s)(4+3s)^2=84+58s.
\]

The exact bounds at both real embeddings give identity and conjugate square-
root uppers 19 and 2, coefficient bound 21, and the 26-element alphabet

```text
0, 1, s, 1+s, 2+s, 1+2s, 3+s, 2+2s, 3+2s, 4+2s,
3+3s, 4+3s, 5+3s, 4+4s, 6+3s, 5+4s, 6+4s, 7+4s,
6+5s, 7+5s, 8+5s, 7+6s, 9+5s, 8+6s, 9+6s, 8+7s.
```

The production left-to-right recursion gives

```text
visited by depth:         1,26,676,3016,7462,12688,117078,493064,425750,156000
Schur-pruned by depth:    0,0,560,2729,6974,8185,98114,476689,419750,117374
exact weight solutions:  2507
nonintegral roots:       1544
integral roots:          963
rejected by old walls:   570
accepted at this shell:  393
```

Among the 570 rejected integral roots, 400 have one explicit lower-wall
blocker and 170 have two.  The manifest freezes the complete root/blocker
witness digest and the blocking-wall occurrence histogram.

## Independent block replay

The hostile replay independently derives the coefficient box, exhausts the
orthogonal `B4`, `F4`, and `A1` component boxes, and performs an exact field-
norm join.  It shares neither the production coordinate order nor its
leading-Schur pruning tree.  Its raw boxes have sizes
`456976,456976,26`; exact bounds retain `488,115,26` vectors, and the keyed
join makes 12,688 outer probes.

It obtains the same 2,507 weights, the same 1,544 nonintegral
reconstructions, the same 570 rejected roots with identical complete blocker
sets, and the same 393 accepted roots.  The reconstructed-root and blocker-
witness digests are

```text
60cf05c9a8618f1154503190b42d23f98b791aad59a53ac35899aaae47699e82
351afd9d5fb1d1f383511d283d7f0552d5044a28d929cb0648f4c1721a474ce7.
```

The manifest stores the full 726-root prefix literally in production order,
as well as all 393 accepted shell roots literally in their production order.
Its deterministic link to the prior 333-wall prefix records the prior
certificate SHA, the exact source reproducer, and the ordered prior-root
digest.  The independent replay verifies equality of the accepted-root set.
Its native lexicographic replay calls the new endpoint wall 443; after root-
vector equality is established, that same endpoint maps to wall 622 in the
certified production order.  Neither Vinberg order is silently sorted.

The replay independently reconstructs the rank-ten inertia witness, full
prefix inertias, Schur endpoint census, and next admissible coordinate.

## Simultaneous admission and prefix inertia

Every accepted root pairs nonpositively with all 333 lower-height walls.  Of
the `binomial(393,2)=77028` distinct same-height pairings, exactly 76,453 are
negative and 575 are zero.  Thus no root in this cohort cuts off another, and
all 393 enter independently of candidate order.

Here “simple” means relative to the exhausted lower-height norm-two prefix and
the complete same-height cohort.  It does not mean that the full chamber wall
list has terminated.

The 726-wall prefix Gram matrix is acute and has exact inertias

```text
identity embedding:   (9,1,716)
conjugate embedding:  (10,0,716).
```

Walls 0 through 9 retain inertias `(9,1,0)` and `(10,0,0)`, hence span the
rank-ten carrier; factorization forces the other 716 directions to be null.

## Local closure of the targeted edge

The diagnostic elliptic edge remains

```text
E = (1,2,4,5,6,7,9,12),   type A2+F4+B2.
```

For 392 new walls, the principal Gram matrix on `E union {r_i}` has identity
inertia `(8,1,0)`, conjugate inertia `(9,0,0)`, and negative identity
determinant.  For new wall 622, the extension instead has inertia `(9,0,0)`
at both embeddings and positive identity determinant.  Thus wall 622 is a
second finite endpoint and the targeted edge's finite endpoint set is

```text
(33,622).
```

This closes only this selected rank-eight edge.  It neither proves that all
other elliptic edges have two endpoints nor audits every affine
zero-neighbourhood.  Those global obligations require a complete simple-wall
system, which is not available here.

## Literal status

The next admissible control coordinate is

\[
x=5+3\sqrt2.
\]

The formal status remains

```text
undecidable_incomplete_simple_roots.
```

There is no proof here of wall completeness or of Vinberg's global compact or
finite-volume criterion, hence no certified new compact or finite-volume
hyperbolic Coxeter polytope.  Infinite volume is not proved either.  Primitive
reflective vectors of norms other than two remain outside this enumeration.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_silver_k3_vinberg_4_plus_3sqrt2.py
PYTHONPATH=src:. python3 scripts/independent_silver_k3_vinberg_4_plus_3sqrt2.py
PYTHONPATH=src:. pytest -q tests/test_silver_k3_vinberg_4_plus_3sqrt2.py
```
