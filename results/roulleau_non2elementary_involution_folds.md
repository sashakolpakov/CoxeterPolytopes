# Two exact Roulleau divisor-chamber folds

This certificate turns two complete divisor-wall configurations in Xavier
Roulleau's [K3 atlas](https://doi.org/10.46298/epiga.2022.6286) into
finite-volume hyperbolic Coxeter polytopes.  It uses
only the Néron--Severi lattice and the complete effective `(-2)`-curve wall
set.  Realization of the wall-graph symmetry by a surface automorphism is not
needed for the polytope proof, but in these two cases it is also known:
Roulleau records `Aut(X)=C2` and marks the hyperelliptic involution as acting
nontrivially on `NS(X)`.  Since the exhaustive weighted wall-graph
automorphism group is itself exactly `C2`, its nontrivial element must be the
folding action certified here.

For the quotient map `pi: X -> X/C2`, finite descent gives
`pi^* Nef(X/C2) = Nef(X) intersect NS(X)^C2` over the reals.  Thus these are
literal quotient-surface nef cones, not merely linear sections with an
unexplained algebraic carrier.  The orbit-sum wall roots are the divisor-side
form of the descended negative curves.

The terminal triangle also supplies the noncommuting folding that the earlier
heuristic left open.  The exchanged curves meet once, so their two simple
reflections generate `A2`; the folded generator is its longest element
`s_1 s_2 s_1`, whose restriction is reflection in their sum, a primitive
`(-2)` vector of divisibility two.  The certificate verifies the resulting
reflection matrix integrally on the whole fixed lattice.  In the terminal
diamond the exchanged curves are disjoint, giving the familiar commuting
product and a primitive `(-4)` folded root.

| source divisor chamber | fixed lattice | polytope | facets | finite vertices | cusps | Vinberg edges |
|---|---:|---:|---:|---:|---:|---:|
| `U+A2(-1)+E8(-1)` (13 curves) | even, det `2`, Smith `1^10 2` | `H^10` | 12 | 18 | 1 | 99 |
| `U+E8(-1)+A3(-1)` (14 curves) | even, det `-4`, Smith `1^10 2^2` | `H^11` | 13 | 26 | 2 | 162 |

For every elliptic edge, Vinberg's criterion finds exactly two finite or ideal
endpoints.  Every full parabolic zero-neighborhood has the required rank, and
there are no bad edges.  Both polytopes are therefore finite-volume and
noncompact.  The source graph automorphism group is exhaustively proved to be
`C2` in each case, the induced action in a divisor basis is integral and
unimodular, the fixed basis is saturated, and every folded root gives an
integral lattice reflection.

The first fold contains eleven roots of type `(norm, divisibility)=(-2,1)`
and one of type `(-2,2)`; its Coxeter diagram has ten `3`-edges and one
parallel pair.  Its Solomon characteristic is
`-1/183936614400`.  The second contains eleven `(-2,1)` roots, one `(-2,2)`
root, and one `(-4,2)` root; its diagram has ten `3`-edges and two `4`-edges.

Independent intrinsic Vinberg searches close the maximal-reflection question,
not just the volume question.  They stop with exactly these walls on
`U+E8(-1)+A1(-1)` at height `1/2` and on `U+D10(-1)` at height `8`.
For each fold the certificate gives an explicit unimodular lattice isometry
that sends every folded wall root to the corresponding intrinsic Vinberg
root.  Thus both chambers have exact index one in the full integral
reflection group.

## Lineage verdict

These are new exact *algebro-geometric constructions in this project*, not
new abstract Coxeter diagrams.  The first normalized diagram agrees exactly
with Scharlau--Kirschmer `List10` entry 3 and two ordinary CoxIter records;
the second agrees with `List11` entry 6 and two ordinary CoxIter records.
Both lie in the already classified `n+2`-facet range of
[Tumarkin](https://arxiv.org/abs/math/0301133), and the pinned
[Scharlau--Kirschmer tables](https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/index.html)
record the corresponding full reflection chambers.
No exact match occurs in the pinned index-two-doubling archive.  Priority for
the specific K3 fixed-section interpretations remains an explicit literature
question; the ledger does not call them new abstract polytopes.

Machine certificate:
`data/enumerations/roulleau_non2elementary_involution_folds.json`.
Rebuild with:

```bash
PYTHONPATH=src:. python3 -m scripts.certify_roulleau_non2elementary_involution_folds --check-archive
```
