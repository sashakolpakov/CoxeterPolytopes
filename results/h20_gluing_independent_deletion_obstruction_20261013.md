# Exact obstruction to deleting walls from the 580-facet H20 patch

**No subset of this fixed dictionary of 580 oriented facet halfspaces
defines a finite-volume Coxeter polytope.** Only two explicit escape
vectors are needed to prove this. This is a finite-dictionary obstruction,
not a nonexistence result for H20 Coxeter polytopes or other constructions.

The dictionary consists of the 580 genuine facets of the
[independently certified six-tile patch](h20_gluing_independent_audit_20261013.md).
Indices below are zero-based in the lexicographically sorted boundary
normal list of the companion boundary census.

## Two walls that cannot be deleted

The common strict timelike center `c` has square `-58` and negative pairing
with all 580 normals. The certificate supplies integer vectors `x_0` and
`x_74` in the ambient facet's orthogonal complement. Their exact data are:

| deleted wall | escape square | pairing with c | pairing with deleted wall | maximum pairing with all other 579 walls |
|---|---:|---:|---:|---:|
| 0 | 22 | -218 | 24 | -2 |
| 74 | 6,718,042,333,490,830 | -1,160,000,134 | 178,461,518 | -2 |

The first witness simplifies to `x_0=c+10*r_0`. The second is an explicit
integer vector in the certificate; its large coordinates do not involve
approximate arithmetic in the verification.

For either deletion, the corresponding spacelike vector satisfies every
remaining inequality strictly. The segment from `c` to this vector
crosses the future light cone while retaining strict inequalities. The
negative pairing with `c` fixes the future orientation of that crossing.
The intersection therefore contains an open set on the ideal sphere and
has infinite hyperbolic volume. Any further wall deletion preserves the
same witness.

Consequently every finite-volume intersection using a subset of these
halfspaces must retain **both** walls 0 and 74.

## The two mandatory walls have a forbidden angle

Their exact Gram is

\[
 \begin{pmatrix}4&-4\\-4&6\end{pmatrix},
 \qquad \cos^2\theta=\frac23.
\]

This is positive definite, and the angle is not `pi/m` for any integer
`m>=2`. In particular it lies strictly between `pi/6` and `pi/4`, so it
could only conceivably be `pi/5`; that angle has squared cosine
`(3+sqrt(5))/8`, different from `2/3`.

There is no facet-redundancy or unrealized-ridge loophole. The certificate
contains an explicit timelike point `p` of square `-364`, with zero
pairing precisely on walls 0 and 74 and strictly negative pairing on
every other one of the 580 walls. Thus these two hyperplanes support true
facets and meet in an actual H18 ridge. Deleting other inequalities
preserves this relative-interior ridge point and its neighborhood.
The same forbidden dihedral angle therefore persists in every enlarged
intersection retaining the pair.

Finite volume forces the pair to be retained, while the Coxeter condition
forbids that pair. This proves the stated obstruction for **all** subsets,
without enumerating the `2^580` choices or relying on a SAT solver.

## Artifacts and replay

The [minimal certificate](../data/enumerations/h20_gluing_independent_deletion_obstruction_20261013.json)
contains the quadratic form, the complete normal dictionary, the common
center, the ridge point, and the two escape vectors. Its
[standalone checker](../scripts/h20_gluing_independent_deletion_obstruction_20261013.py)
uses only Python integers and rational numbers; it needs no LP solver,
NumPy, network connection, or other project verifier.

```sh
python3 scripts/h20_gluing_independent_deletion_obstruction_20261013.py
```

An optional [larger archive](../data/enumerations/h20_gluing_independent_all_deletions_20261013.json)
records exact strict escapes for every one of the 580 single-wall
deletions. It is unnecessary for the two-wall contradiction but can be
replayed with:

```sh
python3 scripts/h20_gluing_independent_deletion_obstruction_20261013.py --all-deletions
```

Before finding this stronger obstruction, the search generated 18
pairwise Gram-compatible subsets, including a symmetric 168-wall subset
and 17 subsets with 360–363 walls. Every one independently received an
exact strict spacelike escape in the
[subset transcript](../data/enumerations/h20_gluing_independent_escape_20261013.json).
Seven exact ideal vertices of the full patch, found by bounded LP
sampling, are recorded in the
[cusp transcript](../data/enumerations/h20_gluing_independent_cusps_20261013.json).
Sampling was only a discovery method; the final obstruction is entirely
the explicit integer certificate above.
