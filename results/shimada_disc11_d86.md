# Exact reduction of the physical (D_{86}) K3 eigensection

Status: complete local half-space reduction and exact traversal across every
internal facet.  This is the physical chamber reached from (D_{68}) across
wall orbit 11, member 1.  It is not yet a global eigen-nef chamber or a
compactness result.

## Exact local result

Shimada's representative (D_{86}) has 35 primitive walls, all in singleton
wall orbits.  After applying the physical transport inherited from the
(D_{68}) crossing and restricting to the fixed ten-dimensional
(\mathbb Q(\sqrt5))-eigenspace, positive proportionality leaves 27 distinct
half-spaces.  Exact face witnesses and positive redundancy identities give

\[
\boxed{
35\text{ ambient walls}
\longrightarrow27\text{ restricted classes}
\longrightarrow11\text{ facets}.}
\]

The eleven facets split as

\[
7\text{ restricted }(-2)\text{-root facets}
\quad+\quad
4\text{ internal Borcherds facets}.
\]

Wall labels below are pairs `(Shimada orbit, member)` in the physical
(D_{86}) state.  The seven root facet classes are

\[
\begin{gathered}
\{(0,0),(13,0)\},\quad
\{(1,0),(4,0)\},\quad
\{(2,0),(16,0)\},\quad
\{(3,0),(14,0)\},\\
\{(6,0),(12,0)\},\quad
\{(7,0),(10,0)\},\quad
\{(8,0),(9,0)\}.
\end{gathered}
\]

The proportionality ratios of the second wall to the first are respectively

\[
\phi-1,\ \phi-1,\ \phi,\ \phi,\ \phi-1,\ \phi,\ \phi-1.
\]

The four internal facets are the singleton classes

\[
\boxed{(20,0),\ (23,0),\ (24,0),\ (26,0).}
\]

For every one of these eleven classes, the certificate gives a point in
(\mathbb Q(\sqrt5)^{10}) that has positive square, pairs to zero with that
class, and pairs strictly positively with all other 33 or 34 ambient walls.
Hence every listed inequality really is a hyperbolic facet; this conclusion
does not come from a floating-point hull computation.

## Exact redundancy

The remaining 16 positive-proportional classes each have an archived positive
(\mathbb Q(\sqrt5))-linear dependence on the eleven facet representatives.
For example, writing (w_i) for the restriction of wall ((i,0)),

\[
w_{11}=(\phi-1)(w_{20}+w_{26})
\]

proves that the root wall ((11,0)) is redundant, while

\[
w_{17}=(2-\phi)(w_{24}+w_{26})
\]

proves redundancy of the class
(\{(17,0),(34,0)\}); the second member restricts to (\phi w_{17}).
All coefficients in all sixteen identities are exactly positive at the
distinguished embedding.

The ambient chamber has 17 `_minustwo` walls.  Exact recovery with the inverse
Néron--Severi Gram matrix proves that all are integral classes of square
(-2).  Under the transported golden structure:

- 16 are integral field-reflection roots;
- their shell distribution is eight of norm (2) and eight of norm
  (4-2\phi);
- their field span has rank 9; and
- wall ((15,0)) is the unique non-field root and is redundant.

Thus the seven root facets are represented by 14 ambient roots.  The other
three ambient root walls, ((5,0),(11,0),(15,0)), do not cut the local
eigensection.

## One star beyond (D_{86})

Every internal facet is a genuine simple crossing: unlike one of the
(D_{68}) facets, none requires a simultaneous multi-chamber passage.
Exact points immediately across the four faces are strictly interior to the
following physical chambers:

| (D_{86}) wall | target representative | target primitive walls | exact crossing scale |
|---:|---:|---:|---:|
| 20 | (D_{55}) | 41 | (1/16) |
| 23 | (D_{30}) | 34 | (1/4) |
| 24 | (D_{27}) | 33 | (1/2) |
| 26 | (D_{68}) | 28 | (1/2) |

For each target point the quadratic norm is exactly positive and every
primitive target-wall pairing is exactly positive.  Wall ((26,0)) restricts
to the negative of the (D_{68}) wall ((11,1)), confirming that this is
precisely the return edge of the physical chamber graph.

The unexplored forward frontier from this branch is therefore

\[
D_{55},\quad D_{27}.
\]

The physical \(D_{30}\) target has since been reduced and proved entrywise
equal to the state reached through the \(D_{181}\) branch; see
`results/shimada_disc11_d30.md`.  The remaining frontier sections must still
be reduced before the global restricted root cell, Coxeter diagram, finite
volume, or compactness can be claimed.

## Reproduction

The external Shimada supplement has SHA-256

```text
fa39e5b59f4a4b2b1b7fe7632e5737d1734ee3eb7e7023ae26bf2b83387afb68
```

Run

```sh
PYTHONPATH=src:. python3 scripts/analyze_shimada_disc11_d86.py /path/to/disc11.txt
SHIMADA_DISC11=/path/to/disc11.txt PYTHONPATH=src:. \
  pytest -q tests/test_shimada_disc11_d86.py
```

The verifier archives all eleven exact face witnesses and all sixteen exact
positive dependences.  It also pins SHA-256 fingerprints of the physical
transport, the 35 transported wall rows, and their exact restricted
functionals, so the certificate is tied to one unambiguous physical copy of
(D_{86}).
