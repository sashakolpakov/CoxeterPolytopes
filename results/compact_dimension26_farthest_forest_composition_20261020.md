# Compact dimension at most26 from actual farthest-vertex projections

20 October 2026. Complete global composition by root. Paper proof;
no scientific, Lean or TeX execution. Acceptance requires the
source-bound audits of every input listed below and a separate
review of this composition. No arbitrary metric-face inheritance
or arithmetic hypothesis is used.

**Theorem, relative to the explicit reviewed component inputs.**
Every compact hyperbolic Coxeter polytope has dimension at most26.

## 1. Exact input boundary

All paths below are under results/. The accepted farthest-vertex
normalization supplies one vector x with ||x||^2<1. Its spatial
components along the irreducible finite vertex diagram are
orthogonal, so their squared norms add. Classical finite spherical
diagram classification and the published ambient compact ceiling29
are external inputs. No classification of whole hyperbolic
polytopes is assumed.

The component inputs are:

| Input | Exact source and SHA256 |
| --- | --- |
| Complete ordinary finite-type rates, small exception accounting, B/D>=9 energy>9/14, and component caps A<=20,B<=22,D<=23 in dimensions<=29 | `compact_forest_ordinary_a_parameter_and_dominant_classical_frontier_20261020.md`, `dd6f8b0de2d7d0bc835be3eb4e320013aa40792407fb2a747f35bcf04385d5ad` |
| Every A_r,r>=7 has actual projection>1/2 | `compact_large_an_actual_half_uniform_transfer_20261020.md`, `7351b471ba3b0edaf0a524ccf082010e458b7fedd8d063491a75632223bd70e7` |
| B_r,r>=21 and D_r,r>=22 have first17-coordinate energy>47/50 | `compact_classical_prefix_energy_and_dimension27_composition_20261020.md`, `5232d17e565d9da173784bb77a43231bb05bb1c5d0042cb685dd7d9d09d42cf3` |
| A20 has actual projection>4/7 | `compact_a20_actual_projection_four_sevenths_20261020.md`, `be7fd5091a83a2fbbdce341bfe7e731b8c96adfbbccf3ed3bb716691c07d2d52` |

For the first row only its component conclusions and type-independent
exception proof are used, not its parameterized A premise or any
then-current interpretation of M. In particular this composition
does NOT require the stronger standalone ordinary A8 theorem.

Every component other than A>=7 or B/D>=9 has squared projection
strictly greater than its rank divided by14, apart from the
following exhaustive alternatives:

* At most one short A1 may be discarded, losing one rank.
* Instead, one isolated rank-two7 defect may pay only>9/100. Then
  there is no short A1 and no second such defect.

The universal small-opposite isolation and native wall-product gap
prove these alternatives with arbitrary spectators. Type J may be
discarded conservatively as one rank, exactly like type I; no
stronger pair-funding claim is needed here. The rates for ordinary
A5/A6, B5--B8, D5--D8 and E6--E8, as well as all lower types,
are included in the pinned first input. All estimates concern the
SAME actual vector x.

## 2. A dominant component exists and is unique

Suppose, towards a contradiction, that27<=n<=29. Call an A component
of rank>=7 or a B/D component of rank>=9 dominant.

If none exists, every component is ordinary except the permitted
small debit. With no rank-two defect, discarding the possible one
short A1 gives

    ||x||^2>(n-1)/14>=26/14>1.

With one rank-two defect it gives

    ||x||^2>9/100+(n-2)/14>=9/100+25/14>1.

Thus a dominant component exists. Every dominant component costs
strictly more than1/2: this is the new A theorem or the already
reviewed B/D rate9/14. Two dominant components would therefore
contradict ||x||^2<1. There is exactly one. All its spectators are
ordinary apart from the same single permitted debit.

Let r be its rank and S=n-r the spectator rank. The three spectator
costs are, when present,

    S/14,  (S-1)/14,  9/100+(S-2)/14.             (1)

When there is no spectator, its energy is0 rather than a strict
positive lower bound. All uses of(1) below have enough positive
ordinary spectator rank, so this empty-sum boundary causes no issue.

## 3. The dominant A case

The dominant A component costs>1/2. If there is no exception,
S>=7 gives total>1/2+7/14=1. With one short A1, S>=8 gives
total>1/2+7/14=1. With a rank-two defect, S>=8 gives
total>1/2+9/100+6/14>1. Consequently S<=7 in every case.

The component cap gives r<=20. Hence n=r+S<=27. Since n>=27,
the only possibility is n27,r20,S7. Now use the stronger actual
A20 rate>4/7. Even after discarding one short A1,

    ||x||^2>4/7+6/14=1.

With a rank-two defect the lower bound is

    ||x||^2>4/7+9/100+5/14=713/700>1.

The fully ordinary spectator case is stronger still. The dominant
A case is impossible.

## 4. The dominant B or D case

The dominant B/D component costs>9/14. The three rows of(1) give
S<=4 with no debit and S<=5 with either debit, because their
contradiction thresholds are respectively S5, S6 and S6. Therefore
r=n-S>=22. The component caps give r<=22 for B and r<=23 for D.
In particular S=n-r>=27-23=4.

The new localized prefix theorem now applies to either type,
giving component energy>47/50. At S>=4 its spectators cost more
than3/14 even after discarding one short A1, or more than
9/100+2/14 with a rank-two defect. Both exceed the remaining
3/50. Thus the total actual squared norm exceeds1, another
contradiction. All possible dominant types are exhausted.

Dimensions27,28,29 are therefore impossible. Together with the
external compact ceiling29 this proves dimension at most26.
This implication concerns arbitrary compact Coxeter polytopes,
not only arithmetic ones; the existing arithmetic ceiling20 is
already stronger and is unchanged.

## 5. Verification and next boundary

This source is the global bookkeeping proof. Its component geometry
belongs to the separately pinned sources and their written audits.
It claims no fresh script execution, Lean compilation, external
refereeing, new example or sharpness. A final conditional Lean
module may assume the geometric rates and complete finite-type
coverage explicitly, then check this integer/rational logic.

This is not a claim that dimension26 occurs. Largest known compact
examples remain in dimension8. Further work should strengthen the
large-A rates and lower-rank B/D prefix estimates enough to reduce
the global ceiling again, rather than expand the parked nine-field
refinement or the already sharp right-angled dimension packages.
