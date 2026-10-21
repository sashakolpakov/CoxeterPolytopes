# Exact Gram obstruction for one Record30 `a=2` orbit-10 labeling

## Result

The locally admissible Record30 `a=2`, positive-`A4` orbit-10 labeling
described below is **not** the Coxeter diagram of a chamber in
\(\mathbb H^9\), for any values of its three dotted entries.

This is a global metric obstruction, not another local-label failure.  The
candidate really does pass the frozen necessary local gates: all 128
rank-nine actual-vertex diagrams are elliptic and all four order-four or
order-five minimal nonfaces are Lanner.  Nevertheless, an exact fixed
principal Gram minor forces rank at least 11.

## Frozen candidate

The inclusion-minimal nonfaces are

```text
(0,11), (1,12), (2,13),
(0,3,4,5), (6,7,8,11),
(3,4,5,9,10), (6,7,8,9,10).
```

The first three pairs are dotted.  In the twice-normalized Gram matrix write
their entries as \(x,y,z<-2\).  Every unlisted meeting pair has Coxeter label
2.  The finite non-right labels are

```text
m=3: (0,1), (0,5), (2,11), (3,5), (3,10),
     (6,8), (6,9), (8,11), (9,10)
m=5: (4,5), (7,8).
```

Thus this is sparse terminal state 0 with the two additional locally allowed
cross-component edges `(0,1)=3` and `(2,11)=3`.

## Decisive exact rank test

Use twice the normalized Gram matrix, so a diagonal entry is 2 and a finite
off-diagonal entry is \(-2\cos(\pi/m)\).  The principal block on facets

```text
0,1,2,3,4,5,6,7,8,9,10
```

contains no dotted entry.  Exact elimination gives

\[
 \det G_{0,\ldots,10}=21-7\sqrt5=7(3-\sqrt5)>0.
\]

The sign is elementary: \(3>\sqrt5\) because \(9>5\).  Consequently this
11-by-11 block is nonsingular, independently of \(x,y,z\), and every full
Gram completion has rank at least 11.  But fourteen wall normals for a
chamber in \(\mathbb H^9\) lie in \(\mathbb R^{9,1}\), so their Gram matrix
has rank at most 10.  The required inertia `(9,1,4)` is therefore impossible.

For an explicit Schur replay, the leading ten-facet anchor has determinant

\[
 52-24\sqrt5\ne0,
\]

and its last Schur residual is

\[
 -\frac{7(9+5\sqrt5)}{44}\ne0.
\]

Their product is \(7(3-\sqrt5)\).

The certificate also derives the determinant independently from the weighted
matching polynomial of the fixed non-right forest.  Its weighted matching
sums in sizes 0 through 5 are

```text
1,
10 + sqrt(5),
35 + 7*sqrt(5),
(99 + 29*sqrt(5))/2,
(51 + 19*sqrt(5))/2,
(7 + 3*sqrt(5))/2.
```

Substitution in

\[
 \sum_M(-1)^{|M|}2^{11-2|M|}\prod_{e\in M}w_e^2
\]

again gives \(21-7\sqrt5\).

## Exact scope

This excludes one exact sparse label state and the exact two-cross-edge
placement `(0,1),(2,11)`, uniformly over all three dotted values.  It does
**not** exclude:

- the other sparse terminal states on this nerve;
- other placements or labels of cross-component edges;
- configurations with three or more cross-component edges;
- the whole orbit-10 nerve; or
- the whole Record30 `a=2` schedule.

Accordingly this is not a new hyperbolic polytope and not yet a schedule-wide
nonexistence theorem.  It also makes no finite-volume noncompact claim.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_h9_14_record30_a2_orbit10_gram.py --verify
PYTHONPATH=src:. python3 -m pytest -q tests/test_h9_14_record30_a2_orbit10_gram.py
```

The frozen machine-readable certificate is
`data/enumerations/h9_14_record30_a2_orbit10_gram_obstruction.json`.
