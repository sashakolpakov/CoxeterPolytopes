# New-example search, 13 October 2026

## Outcome and repository purpose

This repository searches for hyperbolic Coxeter polytopes, with an emphasis
on exact certificates and on advancing the compact frontier beyond dimension
eight. It also studies algebraic realizations of known and newly extracted
noncompact chambers. Those realization results are not automatically new
abstract polytopes.

**This search has not established a new example or a new dimension.**
The computations below are concrete new attempts. A finite prefix that has
infinite volume does not establish nonreflectivity of the complete host:
later mirrors may still close it.

The compact target is H9. Finite-volume examples are known through H19 and
in H21; H20 is another missing dimension. Esselmann's arithmetic obstruction
in dimension 20 concerns **noncompact arithmetic** examples. A noncompact
H20 search must therefore leave that class; compact arithmetic H20 over
nonrational fields is not excluded by that theorem. See the separately
[sourced literature audit](frontier_literature_2026_10_13.md).

## Calibration

The known compact Bugaenko H8 chamber was reproduced and independently
checked: eleven facets, 41 finite vertices, no cusps, and no hyperideal rays.
Its existing geometry and Vinberg tests passed (`9 passed`).

## Quadratic diagonal H9 search

A new exact scan covers 132 coefficient presentations

`L = <-alpha> + <1>^9`, `alpha=a+b*w`,

with `a=-2,...,3`, `b=1,2`, retaining exactly those with `alpha>0` and
`alpha'<0`. The squarefree radicands are

`2,3,5,6,7,10,11,13,17,19,21,29,33`.

Here `w=(1+sqrt(d))/2` when `d=1 mod 4`, and `w=sqrt(d)` otherwise.
The 132 presentations are **not claimed to be pairwise nonisometric**.

For each presentation the search enumerates every norm-one and norm-two
root with Vinberg height `alpha*x^2/q <= 300`, modulo signed coordinate
permutations. It uses the dominant spatial chamber
`y1 >= ... >= y9 >= 0`. The root equation

`sum(y_i^2) = q + alpha*x^2`

bounds every coordinate at both real embeddings. Integer coefficient boxes
are derived from those two bounds; all signs, sums of squares, height
comparisons, and acute-wall decisions use integer quadratic-field
arithmetic. The search inspected 37,223 dominant positive-height roots and
retained 9,926 positive-height walls in total across the presentations.

All 132 displayed prefix cones have an **exact feasible spacelike point**,
so none is a finite-volume example. Each cone has a strict timelike interior
point near the time axis; perturbing a feasible spacelike point toward it
gives a strict spacelike point, proving infinite volume. Simple equal-coordinate witnesses
suffice in 122 cases. For the other ten, bounded numerical linear programs
propose rational points; an integer-only independent verifier checks every
inequality and the strictly positive Lorentzian norm. Numerical success or
failure is never a proof oracle. Other root norm classes have not been
exhausted, and no complete-host nonreflectivity conclusion is drawn.

Artifacts:

- [Search](../scripts/search_quadratic_diagonal_h9_2026_10_13.py)
- [132 archived prefixes](../data/enumerations/quadratic_diagonal_h9_2026_10_13.json)
- [Ten supplementary witnesses](../data/enumerations/quadratic_diagonal_h9_geometry_2026_10_13.json)
- [Standalone exact verifier](../scripts/verify_quadratic_diagonal_h9_2026_10_13.py)
- [Tests](../tests/test_quadratic_diagonal_h9_2026_10_13.py)

Reproduction:

```sh
python3 scripts/search_quadratic_diagonal_h9_2026_10_13.py
python3 scripts/audit_quadratic_diagonal_h9_2026_10_13.py
python3 scripts/verify_quadratic_diagonal_h9_2026_10_13.py
python3 -m pytest -q tests/test_quadratic_diagonal_h9_2026_10_13.py
```

The numerical proposal step requires NumPy/SciPy. The final independent
verifier needs only Python's standard library. The three focused tests
compare arithmetic to the existing exact field implementations, compare
small root searches against a separate Cartesian enumeration, and replay
all 132 witnesses.

## Mixed root norms in the natural golden H9 host

The earlier repository theorem for
`<-2phi> + E8 + <2>` establishes infinite volume for the **norm-two**
reflection chamber. It does not establish that the full integral reflection
group has infinite index. In particular,

`r=(phi-1,0^8,phi)`

has norm four and defines an integral reflection. Its matrix on `(x,u)` is

`[[phi,-1],[phi,-phi]]`,

the same involution that the earlier norm-two analysis used as a chamber
symmetry. Adding its mirror changes the chamber. The new chronological
norm-two/four search follows this distinction explicitly.

The first full mixed-norm prefix, through height 30, has 17 walls and its
exact cone census is 86 finite vertices, zero cusps, and 196 hyperideal
rays. Thus this prefix also fails finite volume. Through height 300 the full mixed-norm search has 53 walls. Its exact
feasible spacelike point has norm `68/9-(14/3)phi > 0`, so that larger prefix
also has infinite volume. The completed continuation through height 1,000
has 215 walls; a new exact escape has norm `233/2-72phi > 0`. The standalone
verifier checked all 215 integral reflections, their chronological heights,
the strict timelike interior point, and every escape inequality. Thus this
deeper prefix also has infinite volume; full-host reflectivity remains open.
The odd variant `E8+<1>` was searched separately
with its full norm alphabet `{1,2}`: through height 300 it has 184 walls,
and its independently checked escape has norm `(61-37phi)/49 > 0`. See the
[dedicated arithmetic report](golden_h9_fresh_search_2026_10_13.md), including
five additional golden path hosts and a known compact H4 calibration.
The new [independent mixed-norm verifier](../scripts/verify_golden_h9_mixed_escape_20261013.py)
uses only standard-library rational arithmetic and checks every integral
reflection, root height, interior inequality, and escape inequality.

## Geometric and nonarithmetic alternatives

Two independent additional lanes were pursued: cap extensions of Bugaenko's
H8 chamber into H9, and exact face extraction from Vinberg's nonarithmetic
H18 mixture. The [cap search](bugaenko_h9_caps_20261013.md) numerically
examined 84,830,744 basis/angle assignments. A retained 34-wall H9 Coxeter
polyhedron passes independent exact angle checks, but a strict spacelike
witness proves its volume is infinite. The numerical dictionary rejection
is expressly not a nonexistence theorem. The latter reconstructs the known 40-wall mixture, then
checks all forty H17 facets: sixteen are Coxeter, but exact Gram matching
identifies all sixteen with known arithmetic chambers. Every seam-crossing
H17 facet fails the angle gate. The
[literature and mixture report](frontier_literature_2026_10_13.md) records
the further face searches and their precisely limited scope: all 12,787
specified seam-crossing H17, H16, H15, and H13 candidates have independently
verified forbidden-angle witnesses. The sixteen positive H17 facets are
exact matches to known chambers; none supplies a new example.


## Validation status

The original Bugaenko calibration has nine passing tests. The existing
71-wall norm-two H9 replay also passes after the backward-compatible
enumerator extension. Ten focused tests cover the new arithmetic, low-rank
root enumeration cross-check, mixed-norm reflections, all 132 diagonal
witnesses, the two H300 mixed-norm witnesses, and rejection of a forged
witness norm. The cap and 12,787 face exclusions additionally have separate
exact verifier runs.

These checks validate the stated bounded results. They do not turn an
incomplete root search or a finite numerical cap dictionary into an
exhaustive nonexistence theorem.
