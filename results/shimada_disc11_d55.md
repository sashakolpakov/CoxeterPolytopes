# Exact reduction of the physical \(D_{55}\) K3 eigensection

Status: complete local half-space reduction and exact traversal across every
internal facet. This is the marked chamber reached along

\[
D_{68}\xrightarrow{(11,1)}D_{86}\xrightarrow{(20,0)}D_{55}.
\]

It is not a claim that the global eigen-nef chamber graph has closed.

## Exact local result

The physical \(D_{55}\) chamber has 41 primitive ambient walls. Restriction
to the fixed ten-dimensional \(\mathbb Q(\sqrt5)\)-eigenspace gives

\[
\boxed{
41\text{ ambient walls}
\longrightarrow33\text{ restricted rays}
\longrightarrow11\text{ facets}.}
\]

All reductions are exact over \(\mathbb Q(\phi)\). Each facet has an archived
timelike relative-interior point, and each of the other 22 restricted rays
has a strictly positive exact conic expression in facet representatives.
Thus the facet assertion is not based on a floating-point hull decision.

The facets split as seven root facets and four internal Borcherds facets.
Using `(orbit, member)` labels in this physical state, the root classes are

\[
\begin{gathered}
\{(0,1),(6,0)\},\quad
\{(1,0),(3,0)\},\quad
\{(1,1),(8,1)\},\\
\{(4,0),(5,0)\},\quad
\{(4,1),(7,1)\},\quad
\{(7,0),(8,0)\},\quad
\{(9,0),(10,0)\}.
\end{gathered}
\]

The four singleton internal facets are

\[
\boxed{(12,1),\ (13,1),\ (14,0),\ (16,1).}
\]

## Divisor/root interpretation

The ambient chamber has 18 `_minustwo` walls. Exact inversion of the
Néron--Severi Gram matrix and the transported golden structure proves:

- 17 walls are integral field-reflection roots;
- their norm shells are eight of norm \(2\) and nine of norm
  \(4-2\phi\);
- their field span has full rank 10; and
- `(3,1)` is the unique non-field root wall, and it is redundant.

Every wall occurring in a root facet is a field root. In particular, this
second local cell reaches full field rank, unlike physical \(D_{86}\), whose
surviving roots span rank 9.

## Exact internal first star

All four internal facets are simple crossings. Exact points immediately on
the other side are strictly timelike and pair positively with every primitive
wall of the target physical chamber.

| \(D_{55}\) wall | target representative | target walls | crossing scale |
|---:|---:|---:|---:|
| `(12,1)` | \(D_{20}\) | 52 | \(1/3\) |
| `(13,1)` | \(D_{14}\) | 40 | \(1/2\) |
| `(14,0)` | \(D_{16}\) | 39 | \(1/3\) |
| `(16,1)` | \(D_{86}\) | 35 | \(1/2\) |

The last line is certified by equality of the complete transported physical
wall sets and of the cumulative row isometries, not merely by the repeated
representative index 86. Moreover, `(16,1)` restricts to the negative of the
entering \(D_{86}\) wall `(20,0)`. The other three targets have pairwise
distinct physical wall sets and do not coincide with the known physical
copies of \(D_{68}\), \(D_{86}\), or \(D_{55}\).

The new forward frontier from this branch is therefore

\[
\boxed{D_{20},\quad D_{14},\quad D_{16}}
\]

as physical marked chambers. Their representative labels alone are not used
as graph keys.

## Reproduction

The external Shimada supplement has SHA-256

```text
fa39e5b59f4a4b2b1b7fe7632e5737d1734ee3eb7e7023ae26bf2b83387afb68
```

Run

```sh
PYTHONPATH=src:. python3 scripts/analyze_shimada_disc11_d55.py /path/to/disc11.txt
SHIMADA_DISC11=/path/to/disc11.txt PYTHONPATH=src:. \
  pytest -q tests/test_shimada_disc11_d55.py
```

The verifier pins the marked state, labelled physical wall rows, unordered
physical wall set, and exact restricted rows by SHA-256 fingerprints.
