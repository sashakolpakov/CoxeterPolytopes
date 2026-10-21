# Two-contact families on the H9 chamber's H3 face

The earlier three-contact bank did not cover continuous walls with only two
prescribed old contacts. This investigation fills that bounded search gap,
solves mutual finite-angle equations, and combines new candidates with the
old bank. It finds **no compact H3 closure**. Its exact mathematical result
is a permanent-cusp obstruction to **23 entire contact recipes**. All other
search results here are numerical diagnostics; no new exact wall,
finite-volume chamber, or construction dimension is claimed.

All files in this package are new. The retained chamber is the independently
reviewed [55-wall chamber](fresh_nonarithmetic_nonright_s6_20261016_cone_review.md).
Its actual face on `S6={0,1,2,3,5,6}` has 13 finite vertices and seven
hyperideal extreme rays. The [projection theorem](fresh_nonarithmetic_h3_escape_projection_20261016_review.md)
already shows that right-S6 walls cannot cut the saved ambient escape.
The purpose here is to close the whole H3 face while retaining the two
nonright-S6 ambient cutters, rather than to claim that local closure would
close the ambient chamber.

## Domain and parameterization

The [old 4,329-root bank](fresh_nonarithmetic_h3_joint_cover_20261016.json)
prescribed three old contacts. In the residual Lorentz four-space this
leaves an affine line, whose norm-two equation ordinarily has isolated
solutions. Two prescribed contacts instead leave a two-dimensional kernel
and a continuous norm conic. This is a genuinely larger domain.

The [new domain source](../scripts/fresh_nonarithmetic_h3_two_contact_20261016.py)
and [saved data](fresh_nonarithmetic_h3_two_contact_20261016.json) use a
future-oriented orthonormal residual basis with metric `diag(-1,1,1,1)`.
Every old normal is first projected into this space; native old squared
norms are not substituted for the projected norms.

The anchor pairs are all 15 pairs occurring among the three non-S6 walls
incident to the seven current exterior rays. Each prescribed label ranges
from 2 through 10. Thus there are 1,215 recipes. For each nonsingular
two-anchor frame, the source solves the affine equations and diagonalizes
the residual two-dimensional quadratic form. It finds 2,268 real conics,
of which 1,144 have numerically admissible open intervals:

| Type | Number | Parameterization |
|---|---:|---|
| Hyperbola branches | 1,039 | `p+a*t+b/t`, `t>0` |
| Circles | 80 | `p+a*(1-t^2)/(1+t^2)+b*2t/(1+t^2)` |
| Null lines | 25 | `p+a*t` |

The circle parameter omits its single projective endpoint; boundary
equalities and other isolated degeneracies are not accepted. No singular
anchor-frame or point-degeneracy record occurred in this bounded pilot.

Remaining nonconstant old products must be strictly dotted; constant right
or finite Coxeter contacts are allowed. The retained future center must
remain strict. Each condition becomes a quadratic polynomial inequality
after multiplication by a positive parameter denominator. Numerical roots
partition the parameter line into open intervals. Before sampling, these
intervals are split at the pairing roots of all seven initial ideal boundary
points, so a narrow initial cutting window is not lost to a coarse global
parameter grid.

This is a numerical interval calculation, not an exact sign certificate.
Other anchor pairs, arbitrary labels, isolated additional finite contacts,
one-contact families and walls with no old finite contact are outside it.
In particular, the finite label range is never used as a universal angle
alphabet.

## Exact removal of the affine line families

The separate [exact gate](fresh_nonarithmetic_h3_two_contact_20261016_affine_gate.md)
uses exact old incidences and three-root Grams. Seven right anchor pairs
with new labels `(3,6)`, `(4,4)`, or `(6,3)`, and the old pi/4 pair `(30,52)`
with new labels `(2,4)` or `(4,2)`, form affine triangles. These are 23
recipes accounting for all 25 numerical null-line branches.

Their positive null vector persists under every acute addition sharing
any strict future interior point. The retained pi/5 pair 0,1 is orthogonal
to each triangle. A finite-volume completion would therefore require a
Euclidean cusp translation lattice preserved by a rotation with irrational
trace, which is impossible. The gate covers every real parameter and every
coefficient field for these recipes. It does not exclude the remaining
1,119 conics. The new searches below remove all 25 branches.

## Coupled finite contacts

The [pair source](../scripts/fresh_nonarithmetic_h3_two_contact_20261016_pairs.py)
and [saved pair search](fresh_nonarithmetic_h3_two_contact_20261016_pairs.json)
choose 160 first-cap samples stratified across all eight observed initial
coverage patterns. Both old anchor labels are at most 6. For each sample
and each eligible second family, every mutual label 2 through 10 is solved
as a quadratic equation in the second parameter, followed by all 55 old
contact gates.

There are 485,280 mutual equations and 55,037 compatible finite-contact
pairs. None covers an entire initial exterior component. The seven exterior
rays form the triangle component `{0,1,3}` and the chain `{2,4,5,6}`.
The saved 24 best bridge pairs cut distinct points only in `{2,4}` or
`{5,6}`. These counts are bounded by the selected first parameters and label
ranges. They do not exclude an unsampled continuous coupled solution.

## Repairing the earlier six-wall trial

The [repair source](../scripts/fresh_nonarithmetic_h3_two_contact_20261016_trial_extension.py)
and [saved repair data](fresh_nonarithmetic_h3_two_contact_20261016_trial_extension.json)
retain the six numerical walls of the earlier
[trial](fresh_nonarithmetic_h3_closure_20261016_trial.json), in addition to
the exact 55. Their norms, mutual compatibility, compatibility with all 55
walls and the old trial rays' inequalities are rechecked numerically.

For every remaining two-contact family, this stage solves finite contact
labels 2 through 10 to each of the six trial walls. It also samples
all-trial-dotted intervals separately in each replacement-ideal coverage
cell. The 60,426 finite equations yield 2,143 new candidate walls. Combining
them with the reconstructed old 4,329-root bank, retaining only cutters of
the trial's exterior samples and deduplicating, gives 4,542 candidates.
Every candidate passes all 61 numerical old-wall gates.

All 10,312,611 candidate pairs are tested. Their incompatibilities are
represented by eight incompatibility cliques and 172 compressed rows.
Each compressed row is itself an incompatibility clique: its two parts
belong to old cliques and every cross pair is incompatible. Conversely every
bad pair is contained in one such row, so this compression does not add a
stronger combinatorial restriction to the computed numerical graph.

A binary program asks for a compatible cover of the sampled ideal boundary.
Each selected collection is clipped against the entire numerical H3 cone;
new exterior rays contribute further ideal boundary samples.

| Iteration | Selected new walls | Finite vertices | Exterior rays | Unresolved null rays |
|---|---:|---:|---:|---:|
| 0 | 5 | 28 | 11 | 3 |
| 1 | 6 | 33 | 7 | 4 |
| 2 | — | Infeasible on 24 boundary samples | — | — |

The last conclusion is the solver's infeasibility result for this finite
floating-point model. It is not a rational dual certificate or a continuous
family exclusion. The near-null rays are explicitly unresolved, not
certified cusps. No cluster is promoted to exact geometry.

## Geometric interpretation and remaining route

Adding compatible walls that cover initial samples can create replacement
exterior regions. Consequently wall counts or initial-ray coverage alone
do not demonstrate progress toward compactness.

A useful necessary condition is intrinsic incidence. For a strict acute
chamber, a positive-definite selected Gram gives an actual face by the
inverse-Stieltjes projection argument. Hence each right-S6 candidate with
positive residual norm is an actual facet of the target H3. If the H3 is
compact, each such polygon facet needs at least three finite neighboring
facets. These intersections must be tested using residual projected Grams;
ambient Coxeter contacts alone do not establish them. This is the same
local geometric criterion explained in the separately developed
[H20 incidence screen](fresh_h20_intrinsic_facet_incidence_20261016.md).

The next broader construction domain is one prescribed old finite contact,
leaving a two-dimensional norm quadric, with the remaining constraints
solved jointly through new-new finite contacts. Another possible domain is
a continuous loop of several two-contact families, without fixing the
first parameter to a sampled value. Neither domain is excluded here.

## Replay and provenance

Generation and fresh saved replay pass for all four sources. The pair
replay initially exposed integer dictionary keys becoming strings in JSON;
explicit string serialization repaired that mismatch without changing any
equations or counts. Source hashes are bound in each JSON. Syntax, direct
and recursive manifests, and local report links are checked separately.

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_nonarithmetic_h3_two_contact_20261016.py --verify
python3 scripts/fresh_nonarithmetic_h3_two_contact_20261016_affine_gate.py --verify
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_nonarithmetic_h3_two_contact_20261016_pairs.py --verify
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_nonarithmetic_h3_two_contact_20261016_trial_extension.py --verify
```
