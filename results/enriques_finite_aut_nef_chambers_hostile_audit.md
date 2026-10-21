# Hostile raw-data audit of the seven Enriques nef chambers

Status: **all independent exact checks pass** (2026-08-08).

## Scope and independence boundary

This is an adversarial replay of the seven finite-automorphism Enriques
full-nef chambers.  The sole chamber input is Shimada's official raw GAP file

```text
data/external/shimada_l10l26note/L10L26compdata2.txt
SHA-256 218c8188e467bde29369cdeda3771d767c3f70105d8fb56def8d698b90e75a6c
```

The hostile script does not import the production Enriques certifier and does
not read its summary JSON or compressed ray archive.  It parses the sixteen
raw `rec` blocks by an independent record-splitting method, chooses its own
lexicographically greedy wall bases, and invokes `cddexec_gmp` directly.  It
uses only the repository's generic exact inertia and literal Vinberg modules.

The Nikulin--Kondō type dictionary is the prescribed published input:

```text
I=12A, II=12B, III=20B, IV=20F,
V=20A, VI=20E, VII=20C (with alternate record 20D).
```

The geometric assertion that these are the complete nodal-curve lists is
also an external Kondō--Nikulin--Martin--Shimada classification input.  This
audit independently checks every lattice, cone, and hyperbolic-volume
consequence of those lists; raw integer arrays alone cannot reprove the
surface classification.

## Independent exact census

For each type the audit verifies

```text
configuration = W * Gram(L10) * W^T,
```

primitiveness and square `-2` of every wall, integrality/involutivity of every
reflection, rank ten, Smith index, normal-Gram inertia `(9,1,f-10)`, and that
every off-diagonal intersection is `0`, `1`, or `2`.  The latter correspond
respectively to Coxeter labels `m=2`, `m=3`, and parallel walls.

| type | raw | facets | wall index | `m2,m3,infinity` | finite | ideal | hyperideal | rank-8 edges | affine audits |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| I | 12A | 12 | 2 | `53,10,3` | 14 | 9 | 0 | 114 | 12 |
| II | 12B | 12 | 2 | `51,15,0` | 38 | 17 | 0 | 270 | 20 |
| III | 20B | 20 | 4 | `138,16,36` | 200 | 44 | 0 | 1,301 | 84 |
| IV | 20F | 20 | 2 | `140,40,10` | 720 | 87 | 0 | 4,405 | 184 |
| V | 20A | 20 | 2 | `135,24,31` | 192 | 29 | 0 | 1,182 | 68 |
| VI | 20E | 20 | 1 | `135,45,10` | 600 | 57 | 0 | 3,420 | 124 |
| VII | 20C | 20 | 2 | `135,30,25` | 450 | 51 | 0 | 2,490 | 92 |

Totals are 2,214 finite vertices, 294 cusps, zero hyperideal rays, 13,182
elliptic rank-eight edges, and 584 affine zero-neighborhood audits.

The seven tuples `(facets,m2,m3,infinity)` are pairwise different.  Hence the
seven normalized Gram matrices cannot be identified by a wall permutation.

## Exact cone and time-sheet audit

For each chamber, ten greedily selected independent walls give coordinates
on the homogeneous cone.  Exact GMP double description enumerates all 2,508
extreme rays.  Every ray is replayed in the ambient `L10` coordinates and
tested against every wall inequality.  Its full active-facet set and exact
Lorentz square are frozen in the hostile JSON.

Every proposed facet has incident-ray span exactly nine, so none of the 12 or
20 walls is redundant.  The sum of all oriented extreme rays is strictly
positive on every wall and has positive square.  Every finite or null ray has
strictly positive Lorentz pairing with this interior vector, proving that all
causal rays lie on one future time sheet.  There are no spacelike/hyperideal
extreme rays and every chamber has at least one null ray.  Thus the exact-cone
test gives finite volume and noncompactness in all seven cases.

## Literal Vinberg audit

Independently of the ray classification, the hostile replay supplies the
complete positive-normal Gram matrices to the literal Vinberg criterion.  It
enumerates every elliptic rank-eight edge, its finite and full-parabolic ideal
endpoint extensions, and every connected affine component's complete zero
neighborhood.

All 13,182 edges have exactly two finite/ideal endpoints.  All 584 affine
zero-neighborhoods are componentwise parabolic of total Coxeter rank eight.
There are zero bad finite-volume edges and zero invalid affine
neighborhoods.  The finite and ideal endpoint counts agree exactly with the
independent cdd census.  Hence all seven decisions are
`finite_volume_noncompact`; finite surface automorphism groups are not being
used as a proxy for the volume test.

## The precise `20C=20D` statement

The raw `20C` and `20D` **configuration matrices are literally equal in the
published ordering**, so they define the same labelled Coxeter Gram and one
Type-VII image.  Their raw 20-by-10 wall-coordinate arrays are **not**
literally equal.  The hostile audit separately verifies for `20D` that its
different wall array reconstructs the common configuration, has primitive
rank-ten walls, has the required `H9` inertia, and spans an index-two
sublattice with Smith diagonal `1^9,2`.

This distinction prevents the shorthand `20C=20D` from being misread as
bytewise equality of the two lattice embeddings.

## Conclusion and claim boundary

No discrepancy with the production seven-image theorem was found.  Given the
published completeness of the curve configurations, the seven raw arrays
independently produce two 12-facet and five 20-facet finite-volume noncompact
full-nef Coxeter polytopes in `H9`.  There is no 14-facet member of this
finite-automorphism Enriques carrier class.  This remains a scoped carrier
statement, not a global nonexistence theorem for `H9/14` polytopes and not a
worldwide novelty claim for the seven known chambers.

## Reproduction and frozen artifacts

```bash
PYTHONPATH=src:. python3 scripts/audit_enriques_finite_aut_nef_chambers_hostile.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_enriques_finite_aut_nef_chambers_hostile_audit.py
```

The fresh regression run, including a complete raw replay, gives
`3 passed in 123.76s`.

- [hostile verifier](../scripts/audit_enriques_finite_aut_nef_chambers_hostile.py)
- [hostile exact certificate](../data/enumerations/enriques_finite_aut_nef_chambers_hostile_audit.json)
- [hostile regression tests](../tests/test_enriques_finite_aut_nef_chambers_hostile_audit.py)

Frozen SHA-256 values:

```text
script       33e97c9ec2f9d82aa4ca66f06bf14879a73a350baa2f8e6f9f61671aaed7fb51
certificate  224ded8c22c223cf379461f54650d2f7fd560292e1de79bbefe62e5d51d55b1b
test         a87050d9f078f91e4b333c76a1557024dd5c6a32986c0753d28dbbd09c02de51
```
