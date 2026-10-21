# Raw hostile audit: original H11 inverse failed as written; explicit repair passes

**Verdict (2026-08-09):**

```text
original primary certificate:  FAIL_AS_WRITTEN
repaired H11 theorem:          PASS
```

The original plane model constructed the two boundaries, ten actual
`(-2)` curves, two crossing `(-1)` curves, and two same-boundary cusp lifts:
`16` of the claimed `26` negative curves.  Its executable then accepted the
hidden-wall theorem by searching for phrases in an upstream JSON.  That is
not a proof that the ten listed but missing cusp-lift classes are effective.

The theorem nevertheless survives.  A helper-free hostile replay constructs
all ten missing curves explicitly over `Q`, proves their absolute
irreducibility and exact cluster multiplicities, and replays the CM and full
cones from the raw vectors.  After this repair, the same-marked-polytope
return is valid.

The independent executable is
[`independent_coble_type_i_two_boundary_plane_inverse_raw_hostile.py`](../scripts/independent_coble_type_i_two_boundary_plane_inverse_raw_hostile.py),
with frozen
[`JSON`](../data/enumerations/coble_type_i_two_boundary_plane_inverse_raw_hostile_audit.json)
and
[`tests`](../tests/test_coble_type_i_two_boundary_plane_inverse_raw_hostile_audit.py).
It imports no project Python helper.  The full replay, including fresh exact
double descriptions and ten `Singular` absolute factorizations, peaked at
`0.086 GiB`.

## 1. Raw-only endpoint theorem

The audit reads only the odd-unimodular Gram and the 26 primitive wall
vectors.  Norms intrinsically recover

```text
2 square -4 boundaries;
10 square -2 walls;
14 square -1 walls.
```

At each contraction step it uses the exact Schur complement, with no integer
truncation.  The exhaustive result is

```text
sequential numerical P2 endpoints:  52
endpoint-degree signatures:          9
line--quintic target endpoints:       4.
```

The four target wall sets are

```text
2,3,5,6,8,10,11,12,13,18,24
2,3,5,6,8,10,11,12,13,19,25
2,4,6,7,9,10,11,12,13,18,24
2,4,6,7,9,10,11,12,13,19,25.
```

For each, the audit derives an exact contraction order and a determinant-one
basis `[H,E1,...,E11]` with Gram `diag(1,-1^11)`.  It then derives, rather
than supplies, the four proximity chains of lengths `3,2,3,3`.  The
degree-one boundary selects `A,B`; the three degree-one roots select the
`BCD`, `AC`-tangent, and `AD`-tangent lines and hence label `C,D`.

Coefficientwise comparison proves, on every branch:

```text
the two raw boundaries = the line and quintic classes;
all ten raw roots      = seven proximity differences plus three root lines;
four contracted -1s   = the four deepest exceptionals;
all 26 plane classes  reconstruct the original raw vectors.
```

After canonical relabeling, all four branches have the same complete
26-class set, not merely the same `2+10+4` visible seed.

## 2. Explicit construction of the ten missing walls

Use the surviving determinantal point `a=b=1` and the canonical affine
normalization

```text
A=(0,0), B=(1,0), C=(0,1), D=(-1,2).
```

For each positive-degree `(-1)` class, the audit independently compiles the
four cluster-jet blocks.  Every resulting linear system has rank one less
than its number of coefficients.  Thus each class has a unique projective
plane curve.  The degree/multiplicity data are

| raw wall on branch 1 | degree | A | B | C | D |
|---:|---:|---:|---:|---:|---:|
| 14 | 2 | `000` | `11` | `111` | `000` |
| 15 | 4 | `000` | `11` | `222` | `111` |
| 16 | 2 | `000` | `11` | `000` | `111` |
| 17 | 4 | `000` | `11` | `111` | `222` |
| 19 | 2 | `000` | `00` | `110` | `111` |
| 20 | 3 | `111` | `00` | `200` | `111` |
| 21 | 5 | `111` | `00` | `311` | `222` |
| 22 | 3 | `111` | `00` | `111` | `200` |
| 23 | 5 | `111` | `00` | `222` | `311` |
| 25 | 2 | `000` | `00` | `111` | `110` |

For example, the four conics are

\[
\begin{aligned}
q_{14}&=x^2+xy-2x+y^2-2y+1,\\
q_{16}&=2x^2+2xy-4x+y^2-4y+2,\\
q_{19}&=3x^2+3xy-2x+y^2-2y+1,\\
q_{25}&=-x^2+xy-2x+y^2-2y+1.
\end{aligned}
\]

The frozen JSON records all ten equations.  Literal blow-up-chart replay
proves that every curve has exactly, not merely at least, its declared
multiplicity at every proper and infinitely-near center.  `Singular`
`absFactorize` returns one absolute component with exponent one for each of
the ten homogeneous curves.

Every strict-transform class satisfies

\[
               E^2=-1,\qquad K.E=-1,\qquad p_a(E)=0.
\]

Absolute irreducibility plus exact transform class and arithmetic genus zero
therefore makes every strict transform a smooth geometrically integral
rational `(-1)` curve.  Together with the four deepest exceptionals, the
plane model now constructs all fourteen `(-1)` walls.

This explicit proof strengthens the parallel H11-specific Riemann--Roch
repair.  The RR proof already establishes effectivity; the interpolation
also gives equations, exact multiplicities, and geometric irreducibility.

## 3. Independent finite-index and cusp gates

From the raw classes alone, the ten actual roots and two crossing curves
reconstruct twelve boundary-orthogonal CM roots.  Their lattice has
determinant `-4`.  Exhausting all half-glues gives exactly one even index-two
overlattice, of determinant `-1`, hence the expected `E10` CM lattice.

A fresh rational double description of the twelve-root chamber gives

```text
extreme rays:                 23
finite / ideal / hyperideal: 14 / 9 / 0
facet face ranks:             9 on all twelve facets.
```

Thus the root reflection chamber has finite covolume; this is an arithmetic
finite-index certificate, not a phrase read from an archive.

The nine primitive CM cusps split exactly as

```text
6 cusps: f-B1/2 and f-B2/2 both integral;
3 cusps: neither lift integral.
```

The twelve integral outputs equal the twelve raw single-boundary wall
vectors coefficient by coefficient.  There are no omitted integral cusp
lifts.

## 4. Exact same-marked-P closure

Kondō's Proposition 3.7 now applies to twelve actual effective irreducible
roots whose reflection chamber has finite covolume.  It excludes hidden
`(-2)` curves and hidden `(-1)` curves meeting both boundaries.

If a remaining hidden `(-1)` curve meets one boundary twice, the standard
primitive-nef-isotropic argument sends it to one of the nine CM ideal rays.
The complete cusp replay therefore gives its class among the twelve raw
lifts.  That class is already represented by one of the explicitly
constructed curves, so it cannot support a distinct irreducible curve.
The strong-Coble negative-curve trichotomy leaves no other possibility.

Hence both implications are proved:

```text
each of the 26 raw wall classes is an actual negative curve;
each actual irreducible negative curve has one of those 26 classes.
```

Therefore the projectivized nef cone is literally the input marked chamber.
A second fresh double description of the full raw cone gives

```text
extreme rays:                   550
finite / ideal / hyperideal:    508 / 42 / 0
face rank on every one of 26 facets: 11.
```

So the polytope is finite-volume and noncompact, independently of the
same-P argument.

## 5. Scope and lineage

The repaired result is a project-new raw-wall derivation of a plane
line--quintic model and an exact inverse for the already frozen `H11/26`
chamber.  It is not a new Coble surface, a new abstract Coxeter polytope, a
compact example, or a class-wide L4 theorem.  Worldwide priority for this
specific plane model remains open.

Primary geometric input: Shigeyuki Kondō,
[*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591),
especially Proposition 3.7 and the primitive-isotropic/genus-one results.

## 6. Bounded replay

```bash
PYTHONPATH=src:. python3 scripts/run_with_rss_limit.py --memory-gib 0.5 --timeout 120 -- \
  python3 scripts/independent_coble_type_i_two_boundary_plane_inverse_raw_hostile.py \
  --check-output data/enumerations/coble_type_i_two_boundary_plane_inverse_raw_hostile_audit.json

PYTHONPATH=src:. python3 scripts/run_with_rss_limit.py --memory-gib 0.5 --timeout 120 -- \
  python3 -m pytest -q \
  tests/test_coble_type_i_two_boundary_plane_inverse_raw_hostile_audit.py
```

