# A direct H9 cap-construction search

**Outcome (2026-10-13): no compact or finite-volume example found.** An
explicit 34-wall candidate passes the exact Coxeter-angle gate, but a
separate exact verifier proves that its volume is infinite. The numerical
failure of the larger finite search dictionary is not a nonexistence
theorem for H9 or for this construction method.

This experiment tries to increase dimension directly from the calibrated
compact eleven-wall Bugaenko H8 polytope. It differs from continuing the
infinite-volume norm-two chamber of the fixed golden H9 lattice: individual new wall normals
may involve additional positive square roots, with no fixed integral H9
lattice imposed.

## Construction and search

Let the eleven archived Bugaenko roots be `r_i`, with square two, in a
Lorentz space of signature `(8,1)`. Add one positive coordinate. Keep the
eleven side walls `(r_i,0)` and the base wall `(0,-sqrt(2))`. A candidate
cap has normal

\[
 (v,\sqrt T),\qquad T=2-(v,v),\quad(v,v)<0.
\]

It has square two and is disjoint from the base. Its side pairings are
exactly `(v,r_i)`. Nine prescribed pairings determine `v`; the remaining
two side pairings must also pass the angle gate. Caps are selected only
when they have mutually allowed Coxeter pairings.

The first numerical dictionary tries labels `2,3,5` at nine independent
side walls, allowing either a finite angle from the same alphabet or an
ultraparallel relation at the other two walls. Of the 55 choices of nine
side roots, 43 are independent. The search evaluates
`43 * 3^9 = 846,369` assignments and retains 130,584 distinct numerical cap
directions after the timelike and side-angle gates.

A compatibility-constrained covering search initially covers 1,041 horizon
samples, including the 41 original vertical vertex directions. Full
numerical cone enumeration then supplies additional escape directions:

| horizon constraints | candidate walls | cone rays | hyperideal rays |
|---:|---:|---:|---:|
| 1,041 | 30 | 1,975 | 1,464 |
| 1,081 | 34 | 2,934 | 2,048 |

Adding the next 40 directions makes this finite numerical dictionary
incompatible. The search then expands the nine prescribed labels to
`2,3,4,5,6`: `43 * 5^9 = 83,984,375` assignments are evaluated. This
expansion retains only caps with projected norm `-100 < (v,v) < 0` that
cut at least one of the 41 original directions or the 80 new escape
directions. It contributes 6,366 retained entries before comparison with
the earlier stock. Their addition does not produce a compatible cover.

These census and cone counts are **numerical exploratory results**:
floating-point angle tolerances, eight-decimal deduplication, Qhull, and
the selected finite dictionaries are not exact completeness certificates.
Other cap parameters, more disjoint side walls, other angle labels,
different base facets, and different H9 constructions remain open.

## Exact rejection of the 34-wall candidate

The retained candidate has 22 caps. Its archived recipes use only exact
pairings `0,-1,-phi`, where `phi^2=phi+1`. The independent verifier
reconstructs every projected vector over `Q(phi)` directly from the
calibrated Bugaenko Gram matrix. It checks the Lorentz signature, every
side-cap angle, and all 231 cap-cap pairs. The latter have labels

| relation | count |
|---|---:|
| 2 | 66 |
| 3 | 19 |
| 5 | 15 |
| ultraparallel | 131 |

The two positive square roots in a cap pairing require no numerical
arithmetic. If the pairing is `b + sqrt(T_i*T_j)`, comparison with a
target `c` is reduced to the sign of `c-b` and the exact comparison of
`(c-b)^2` with `T_i*T_j`.

The certificate supplies rational coefficients, in the basis of the
first nine Bugaenko roots plus the new coordinate, for both a strict
timelike point and a strict spacelike point. The latter has exact square

\[
 \frac{4389895635742039857}{100000000000000000}
 -\frac{13552384375827957}{500000000000000}\phi>0.
\]

Every one of its 34 inequalities is strictly satisfied. For a cap, this
is checked by writing its pairing as `a+z*sqrt(T)`, checking `a<0`, and
proving `a^2-z^2*T>0`, all over `Q(phi)`. The timelike witness is strict
as well. The segment between these two points therefore meets the light
cone strictly inside every halfspace. An open patch of the ideal boundary
survives, proving infinite hyperbolic volume.

Thus the 34-wall system is a valid infinite-volume Coxeter polyhedron,
not the finite-volume or compact new example sought by the project.

## Reproduction

The exact verification needs only the repository's usual dependencies:

```sh
PYTHONPATH=src python3 scripts/verify_bugaenko_h9_caps_20261013.py
```

It returns `PASS_EXACT_INFINITE_VOLUME`. The numerical exploration needs
NumPy and SciPy:

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/explore_bugaenko_h9_caps_20261013.py --expand
```

Artifacts:

- [Exact recipes and witnesses](../data/enumerations/bugaenko_h9_caps_20261013.json)
- [Independent exact audit](../data/enumerations/bugaenko_h9_caps_20261013_exact_audit.json)
- [Numerical search transcript](../data/enumerations/bugaenko_h9_caps_20261013_search.json)
- [Exact verifier](../scripts/verify_bugaenko_h9_caps_20261013.py)
- [Exploratory search](../scripts/explore_bugaenko_h9_caps_20261013.py)
