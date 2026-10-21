# An infinite path-times-cube family in `H12`

## Result

The full-Picard Jacobian--Kummer `H12/172` chamber contains eight exact
configurations of five doubling walls whose reflection group factors as

\[
  D_\infty\times C_2^3.
\]

For the canonical configuration, the parallel pair is W68/W63 and the
orthogonal cube is W31/W58/W59.  Every finite path of `N>=1` chambers in the
first factor, multiplied by the complete three-cube, is a convex
finite-volume noncompact Coxeter polytope made from `8N` seed chambers.  The
exact all-parameter theorem gives

\[
  f(N)=886N+55
\]

for every integer `N>=1`.  Both path phases and all eight listed
configurations give one lattice-isometry class for fixed `N`.  Thus the
canonical sequence begins

```text
N:       1     2     3     4     5     6     7     8
facets: 941  1827  2713  3599  4485  5371  6257  7143.
```

The existence, finite volume, noncompactness, and volume multiplier `8N` are
all-parameter consequences of iterated even-angle doubling.  A separate
no-import hostile audit proves the displayed facet formula for all `N` by an
exhaustive symbolic solve over every integral translation exponent; the
production archive redundantly reconstructs both phases of all eight
markings through `N=8`.

The production builder is
[`certify_jacobian_kummer_h12_path_cube_family.py`](../scripts/certify_jacobian_kummer_h12_path_cube_family.py),
its archive is
[`jacobian_kummer_h12_path_cube_family.json`](../data/enumerations/jacobian_kummer_h12_path_cube_family.json),
and the regression gate is
[`test_jacobian_kummer_h12_path_cube_family.py`](../tests/test_jacobian_kummer_h12_path_cube_family.py).

## The eight configurations

The one-based path pair and cube walls are:

```text
(76,72) x (36,39,48)   (75,71) x (36,42,51)
(74,70) x (39,42,56)   (73,69) x (48,51,56)
(68,63) x (31,58,59)   (67,64) x (32,58,60)
(66,65) x (33,59,60)   (62,61) x (31,32,33).
```

For every row, the two path roots have Gram matrix

\[
 \begin{pmatrix}-4&4\\4&-4\end{pmatrix},
\]

the cube Gram is `-4*I3`, and every path--cube pairing is zero.  All five
walls belong to the complete 124-wall doubling census.  Their integral
reflections therefore give a faithful direct product, with no numerical
angle inference.

The independent centralizer computation proves that the eight marked
configurations form one orbit, with stabilizer 96, and that 48 stabilizer
elements swap the path endpoints.  Hence all eight markings and both phases
are isometric for fixed `N`.  No exhaustiveness claim is made for other
path-cube configurations in the 124-wall stock.

## Exact boundary reconstruction

For the canonical phase, the builder forms all `8N` integral lattice
isometries, removes precisely those selected mirrors shared by two retained
copies, and identifies equal unoriented primitive boundary roots.  For every
`1<=N<=8` it records the sorted root digest, strict interior square, norm and
divisibility histograms, and a rank-13 certificate modulo `1000000007`.
Every boundary root is the image of a primitive seed root under an integral
involution product.  Hence it is primitive, and its reflection is the
integral conjugate of the corresponding seed reflection.

The first two rows agree root-for-root with the independently frozen cube
archive:

- `N=1` is the W31/W58/W59 cube with 941 facets;
- `N=2` is the full W31/W58/W59/W68 cube with 1,827 facets.

The primitive norm histograms continue linearly; for example, `N=8` has
`120*(-2)+6039*(-4)+984*(-8)` boundary roots.

## Independent all-integer audit

The
[`hostile path-times-cube audit`](jacobian_kummer_h12_path_cube_hostile_audit.md)
imports no project module and reconstructs the lattice and 172 walls from the
pinned Dolgachev--Shimada data.  After removing the two path endpoints, it
splits the 939 cube-side walls as `17+36+36+850`.  It then solves 3,526,884
ordered signed side-wall equations for **all** integer translation exponents.
The only finite positive collisions are 922 at `T^0`, 36 at `T^0*S_right`,
and 36 at `T^-1*S_right`; there are 34 all-integer identities and zero
negative, cross-facet, or endpoint--side collisions.  Consequently

```text
17 + 36(N-floor(N/2)) + 36(N-floor((N-1)/2)) + 850N + 2
  = 886N+55
```

for every `N>=1`.  Deleting either endpoint from the 941-wall cube boundary
has Smith diagonal `1^13`, proving full Picard span for every descendant.
The same audit enumerates the order-23,040 published chamber group and its
order-768 centralizer, proving the conjugacy statement above.

## Volume and algebraic carrier

The seed already has independent literal Vinberg and exact-cone
finite-volume certificates.  Every generator is a doubling wall, and the
path and cube factors commute.  The union is therefore an iterated convex
Coxeter doubling with disjoint seed interiors.  Its volume is exactly `8N`
times the seed volume.  Since its reflection group has finite index `8N` in
the nonuniform seed group, every member remains nonuniform and its polytope
is noncompact.

All walls remain integral divisor classes in the same rank-13 full Picard
lattice

\[
 U(2)\oplus D_5(-1)^2\oplus\langle-4\rangle.
\]

Thus the entire construction stays on the seed's seven-dimensional K3
moduli family as a weak divisor-reflection carrier.  Doubling crosses seed
nef walls; this is not a full-nef, effective-curve, or surface-automorphism
statement.

## Novelty boundary

The family and its counts are new to the project ledger.  Coxeter doubling
and finite-path gluing are classical; the positive project content is this
exact five-wall product inside the Jacobian--Kummer carrier and its certified
large boundary systems.  Worldwide priority remains open pending a wider
literature audit.  Absence by facet count from the pinned ordinary and
index-two CoxIter snapshots is bounded evidence only.

## Replay

```bash
PYTHONPATH=.:src python3 -u scripts/certify_jacobian_kummer_h12_path_cube_family.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_jacobian_kummer_h12_path_cube_family.py
python3 scripts/audit_jacobian_kummer_h12_path_cube_hostile.py --check-archive
pytest -q tests/test_jacobian_kummer_h12_path_cube_hostile_audit.py
```
