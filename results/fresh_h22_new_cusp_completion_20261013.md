# H22: completing a full orbit of new cusp links

Three exact placements of a new wall orbit complete the **241,920 specified
deficient cusp links** in the preceding H22 construction. Each resulting chamber
has **5,397,548 actual facets**, and each still has **infinite volume**, proved by
an exact strict spacelike escape. The closest placement cuts the escape that
survived the first two. This constructs compatible cusp repairs and improves the
boundary locally, but does not establish a new finite-volume dimension.

The construction fixes the preceding new angle label at `k=6`; it works for
every old center label `m>=3`. The three placements are separate chambers, each
retaining the same 2,736,428-facet core. Their new walls are replacements for one
another; no compatibility assertion is made about their union.

## A complete local repair

The [preceding face-orbit construction](fresh_h22_gap_face_family_20261013.md)
has new normals `(p,sqrt11)` of norm `2/3`. Each has exactly one parallel partner
`(p',sqrt11)`. Their sum is a future null vector

```text
delta = (p+p', 2sqrt11),
p² = p'² = -31/3,    (p,p') = -35/3,
(p+p')² = -44,       (h,p+p') = -72,    h² = -14.
```

At a representative cusp, the complete incident old wall set splits into seven
orthogonal spherical components:

| Old wall indices | Type |
|---|---|
| 1,21,28,39,46 | A5 |
| 2,3,23,50 | A4 |
| 9,22,45 | A3 |
| 10,34,48 | A3 |
| 32 | A1, norm 2 |
| 63 | A1, norm 2 |
| 198 | A1, norm 4 |

For each component let `a_i` be the negative sum of its simple normals. These
are negative highest roots; the nontrivial components are simply laced chains.
Their norms are `2,2,2,2,2,2,4`. Two additional rational old-space vectors `v,w`
are perpendicular to these eighteen normals and to both `p,p'`. Orthogonalizing
an exact nullspace basis gives

```text
v² = 279/14,   w² = 55/62496,   (v,w) = 0.
```

Use eleven base vectors, in the fixed order

```text
a_0,a_1,a_2,a_3,a_4,a_5,a_6, v,-v,w,-w.
```

For a positive shift `s_i`, add the normal

```text
r_i = (a_i + s_i*(p+p'), 2*s_i*sqrt11),
```

where the same formula uses `v,-v,w,-w` for the last four bases. Its norm is
unchanged, and it is perpendicular to `delta`. Reflect the added coordinate for
the negative-height cusp. At either sign the eleven caps, the eighteen old
walls, and the original parallel pair give exactly 31 incident walls with link

```text
affine A5 + affine A4 + 2 affine A3 + 6 affine A1.
```

The component ranks sum to 21, the horosphere dimension. Each old A component
and its new highest-root cap has a positive all-one affine relation, with sum
`s_i*delta`. The original pair sums to `delta`; the two residual pairs sum to
`(s_7+s_8)*delta` and `(s_9+s_10)*delta`. Together with a strict interior point,
these positive relations prove a bounded Euclidean link. Rank alone would not
give this conclusion.

## Transporting the repair around the whole orbit

There are 120,960 distinct old projections `p+p'`, and 241,920 signed cusps.
The [frame array](../data/enumerations/fresh_h22_new_cusp_completion_20261013_frames.npz)
stores their integer numerators together with both transported tangents. The
exact checker verifies closure under all six old chamber generators,
connectivity from the representative, and equivariance of both tangent frames.
In particular, different generator paths produce the same transported vectors.

The [highest-root array](../data/enumerations/fresh_h22_new_cusp_completion_20261013_highest_roots.npz)
has seven rows per cusp. Its entries are independently reconstructed from the
complete incident subsets of the original 210 walls at every cusp, including
their chain diagrams and norms. Thus there are 1,330,560 new positive-height
walls and 2,661,120 signed new walls. Adding the core gives 5,397,548 facets.

All three placements pass every pairing gate. Eleven representatives are paired
against each entire new orbit; symmetry transfers these comparisons to all
pairs. Pairs at distinct cusps are strictly dotted. At a single cusp the finite
and parallel relations are the affine diagrams above. Opposite-height new pairs
are dotted: the old projections are future timelike, and reversing one height
decreases their pairing. Complete old-wall comparisons and rigorous bounds
against all earlier cap families cover the retained core. Native integer
calculations have explicit overflow bounds.

The old center is strict for every wall. Nonpositive pairings ensure its
orthogonal projection onto any one wall is timelike and strict for every other
wall, proving each inequality defines an actual facet. Strict nonlocal dotted
pairs and the distinct local normals also exclude duplicated walls. The new
old-space projections are timelike, so all new walls are strict at the old
inherited cusps; their previously completed links are preserved.

## Three placements and the remaining escapes

The eleven shifts, ordered as above, are:

| Placement | Shifts | Exact squared norm of a strict escape |
|---|---|---|
| Uniform | `10,10,10,10,10,10,10,10,10,10,10` | `366040322341/10^12` |
| Closer | `3,3,4,4,4,4,4,10,5,1/20,1/20` | `366040322341/10^12` |
| Finite-angle placement | `3,3,4,4,1,1,2,10,5,1/20,1/20` | `3511216341/10^10` |

The [uniform certificate](fresh_h22_new_cusp_completion_20261013_verify.json)
also verifies that the new walls cut the preceding cusp-tangent escape. A
different rational point, given in the
[uniform escape certificate](fresh_h22_new_cusp_completion_20261013_escape.json),
is strict on all 5,397,548 facets. The same point survives the
[closer placement](fresh_h22_new_cusp_completion_20261013_tight.json), as checked
by its [separate escape replay](fresh_h22_new_cusp_completion_20261013_tight_escape.json).

The [finite-angle placement](fresh_h22_new_cusp_completion_20261013_lower.json)
moves the three singleton caps to additional permitted contacts with old walls.
Each moved norm-two representative has eighteen old right angles, one old
`pi/3` relation and two old parallel relations; the moved norm-four cap has
seventeen old right angles, three old `pi/4` relations and one old parallel
relation. Every remaining old relation is dotted.

This placement cuts the previous exact remote escape. For instance, the cap
at highest-root entry 613840, cusp index 87691, has pairing

```text
(-1399991384400 + 425936800800*sqrt11) / 260400000000 > 0
```

with that point. The positivity is checked by an integer-square comparison.
Nevertheless a new rational strict spacelike point exists, of square
`3511216341/10^10`. Its coordinates, all family maxima, and the target cut are
saved in the [lower-placement escape certificate](fresh_h22_new_cusp_completion_20261013_lower_escape.json).
Every radical sign decision is exact. Both height signs and every old center
label `m>=3` are covered. A segment from the common strict timelike center to
such a spacelike point exposes an open part of the ideal boundary, proving
infinite volume in each case.

The bounded-link conclusion is limited to the 241,920 specified cusps and the
previously inherited links. A follow-up
[exact incidence certificate](fresh_h22_lower_new_cusp_20261013.json) diagnoses
two additional null rays: each norm-two singleton cap is also parallel to old
wall 26. At each resulting cusp, exactly seventeen original walls and that one
cap are incident; every other wall in the full chamber is strict. The link has
finite components `A5+A4+2A3+A1` and one affine A1, quotient rank 17 and four
explicit independent lineality directions. Thus new unbounded links remain.
The [source review and fresh replay](fresh_h22_lower_new_cusp_20261013_review.md)
pass. There is no exhaustive census of all cusps or vertices here, and no
exclusion of arbitrary other compatible additions. The following separate
incidence calculation excludes a specified larger family with one cap fixed.

## A conditional obstruction for every real shift in the same families

A subsequent [complete incidence-family census](fresh_h22_lower_incident_family_20261013.md)
fixes the lowered cap `c_4=delta_0-w_32` and all original 210 facets. It then
allows every real shift in all eleven transported affine-ray types, both
height signs, and both normal orientations. This is a larger search than the
three placements above, although it keeps `c_4` fixed.

At the exposed cusp `delta_4=c_4+w_26`, a genuinely new acute-compatible
incident wall must be orthogonal to both `c_4` and `w_26`. These give two linear
equations in its shift. No orbit direction makes both shift coefficients zero,
so all real possibilities reduce to a finite list of rational solutions. The
strict-center condition loses no admissible incident wall: the center is a
positive multiple of the sum of all 210 old normals, and equality would force
the old projection to vanish, contradicting orthogonality to `c_4`.

The complete signed census has 1,126,349 incidence solutions. It excludes 506
zero-center solutions and separately accounts for 184,320 literal copies of
old incident facets. Such copies add no link constraint; copies of `c_4` or
`w_26` are likewise ignored before imposing the genuinely-new-wall equations.
Only five literal candidates, representing three new normals, survive the
seventeen incident-old-wall tests. All three have a positive pairing with a
retained original facet:

| Base type | Required shift | Positive old-wall pairing |
|---|---|---|
| Norm-four highest root | `1/5` | `1/5` with wall 8 |
| First tangent, norm `279/14` | `3/10` | `333/140` with wall 2 |
| Second tangent, norm `55/62496` | `-1/372` | `11/7812` with wall 0 |

The finite-angle gate has no arbitrary label cutoff. Every solved normal has
rational pairings and norm in the rational coordinate form `diag(Q,11)`, so
its normalized squared pairing with an old wall is rational. The only finite
labels are therefore 2,3,4,6. Dotted relations cannot occur between positive
normals in the common null-orthogonal space; its Euclidean quotient gives the
Cauchy bound. Parallel pairs are included explicitly.

Thus no genuinely new acute-compatible wall from these families can be incident
at this cusp. Every finite addition leaves its deficient link, and hence
infinite volume, while all original facets and `c_4` remain facets. This does
not exclude replacing or removing those facets, changing `c_4`'s shift, or
using different normal families. The precise retained-facet condition matters:
the proof is not an assertion about arbitrary inequalities that could make
old walls redundant.

The [independent source review and fresh replay](fresh_h22_lower_incident_family_20261013_review.md)
pass, including a separate Fraction check of the three final rejection
witnesses. The root task also read the final source and freshly replayed the
complete saved census successfully.

## Verification and numerical search scope

The [independent review](fresh_h22_new_cusp_completion_20261013_review.md)
replays all three exact constructions and all three escape checks. Separate
Fraction arithmetic reconstructs the 31 incident normals, their 465 pairs,
positive affine relations, quotient rank, all 2,310 representative/old-wall
relations, and the lower target cutter. It also checks the new escape's norm,
orientation and old-wall signs independently. The large distant-pair census is
a fresh replay of the primary implementation, not a second implementation.

The holonomy and orbit scouts build proposal arrays; the exact checkers verify
their content and do not import numerical search code. Linear programming
proposes remote escapes, which are perturbed and rounded before exact replay.
The [lower-branch scout](../scripts/fresh_h22_new_cusp_completion_20261013_lower_branch_scout.py)
tests ten selected finite-angle anchor branches. Its successful rational
singleton choices are the ones promoted above. That scout is not an exhaustive
label or placement classification.

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h22_new_cusp_completion_20261013_verify.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h22_new_cusp_completion_20261013_escape.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h22_new_cusp_completion_20261013_tight.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h22_new_cusp_completion_20261013_tight_escape.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h22_new_cusp_completion_20261013_tight.py --lower --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h22_new_cusp_completion_20261013_lower_escape.py --verify
python scripts/fresh_h22_new_cusp_completion_20261013_local_review.py --lower --verify
python scripts/fresh_h22_new_cusp_completion_20261013_lower_escape_review.py --verify
```

The preceding exact compatibility and signature certificates are explicit
hash-anchored dependencies. The complete cusp-orbit repairs are established;
finite volume remains the unresolved constructive requirement.
