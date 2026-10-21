# Shimada `X_1`: complete source-stabilizer no-go above `H9`

## Result

Shimada's theorem-complete induced chamber `D^(0)` on the singular complex
K3 surface `X_1` produces no hyperbolic Coxeter polytope above `H9` from a
fixed section by a subgroup of the genuine ample stabilizer

\[
\operatorname{Aut}(X_1,a_1)\cong \operatorname{PGL}_2(\mathbb F_9).
\]

This is an exact no-go in a specified 720-element source-group denominator.
It is not a no-go for the full abstract weighted-wall graph group, if that
group is larger, and it makes no inference about other chambers or K3
surfaces.

## Complete algebraic stock

The primary source is Shimada, [*The automorphism groups of certain singular
K3 surfaces and an Enriques surface*](https://arxiv.org/abs/1412.6904).
Its theorem `thm:D0` proves that `D^(0)` is an induced chamber in the nef cone
and that its automorphism translates tessellate the nef cone.  Theorem
`thm:orbits` and `table:main` give the complete wall set.

The pinned supplement contains:

- the rank-20 Néron--Severi Gram matrix of determinant `-24` and signature
  `(1,19)`;
- the ample class `a_1`, of square 30;
- 1,094 primitive wall vectors, in fourteen orbits of sizes
  `45,45,10,30,72,60,12,40,120,120,120,180,120,120`; and
- all 720 integral matrices of `Aut(X_1,a_1)=PGL_2(F_9)`.

Every source hash, wall norm, ample pairing, Gram isometry, and ample-fixing
condition is replayed by the stock ingester.

## Exact high-subgroup denominator

The fixed-rank histogram of all 720 matrices is

| fixed rank | 2 | 4 | 6 | 8 | 10 | 12 | 20 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| elements | 144 | 324 | 90 | 80 | 36 | 45 | 1 |

Thus the 36 rational involutions give only `H9` and fall outside the strict
above-`H9` queue.  The 45 rank-12 symplectic involutions form one conjugacy
class: an exact commutator census gives centralizer order 16, hence class
size `720/16=45`.

The exact high-pair closure finds 30 `V4` subgroups.  Every one has common
fixed rank 8.  Consequently a proper subgroup with common fixed rank at least
11 is conjugate to the unique symplectic `C2` class.  This is a complete
subgroup argument inside the stated genuine surface group, not an
element-only heuristic.

## Complete restriction and Coxeter obstruction

For a representative symplectic involution, all 1,094 inequalities form 594
wall orbits.  Orbit sums give 580 distinct oriented rays:

| restricted square | negative | isotropic | positive |
|---|---:|---:|---:|
| count | 162 | 28 | 390 |

Every negative restriction has the elementary exact relative-interior
witness obtained by orthogonally projecting `a_1` to its mirror.  Each
witness is timelike, lies on its target mirror, and is strictly positive on
all other 579 restrictions.  Hence all 162 negative restrictions are genuine
irredundant hyperbolic facets.

The exact complete pair gate checks all `C(162,2)=13,041` pairs.  It finds
915 forbidden Coxeter products.  The first independently replayed pair has
facet squares `-4,-6`, inner product `4`, and normalized square

\[
\frac{4^2}{(-4)(-6)}=\frac23,
\]

which is neither a finite Coxeter cosine square nor parallel/ultraparallel.
Both facets have separate strict timelike witnesses.  One such pair already
proves that the fixed cone is not a hyperbolic Coxeter polytope; the 915 count
is the exhaustive diagnostic.

No Vinberg endpoint or cone-volume computation was run, because volume is
irrelevant after the Coxeter-angle prerequisite fails.  In particular this
is a non-Coxeter theorem, not a claim of infinite volume.

## Scope and next stock

The frozen denominator is exactly the theorem-supplied genuine surface group.
A full automorphism computation of the dense 1,094-node oriented weighted
wall graph was deferred while the Record-31 compact search had memory
priority.  It must be performed before the entire abstract chamber-symmetry
denominator can be called closed.

The strongest adjacent unmined source is Shimada's `X_0` chamber from the
same supplement: 3,278 complete walls and a 1,440-element genuine ample
stabilizer.  Its exact element-rank pre-screen has 45 rank-12 nonidentity
elements and no rank-11 class, so its first bounded gate is again a single
`H11` class before any weighted-graph enlargement.

## Reproduction

```bash
PYTHONPATH=.:src python3 scripts/ingest_shimada_aut3singk3_x1_1094_stock.py --check-archive
PYTHONPATH=.:src python3 scripts/certify_shimada_aut3singk3_x1_symplectic_h11_no_go.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_shimada_aut3singk3_x1_symplectic_h11_no_go.py
```

Frozen artifacts:

- `data/enumerations/shimada_aut3singk3_x1_1094_wall_stock.json`
- `data/enumerations/shimada_aut3singk3_x1_symplectic_h11_no_go.json`
- `scripts/ingest_shimada_aut3singk3_x1_1094_stock.py`
- `scripts/certify_shimada_aut3singk3_x1_symplectic_h11_no_go.py`
- `tests/test_shimada_aut3singk3_x1_symplectic_h11_no_go.py`

