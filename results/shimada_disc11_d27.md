# Exact reduction of the physical \(D_{27}\) K3 eigensection

Status: complete local half-space reduction and exact traversal across every
internal facet. The marked chamber is reached by

\[
D_{68}\xrightarrow{(11,1)}D_{86}\xrightarrow{(24,0)}D_{27}.
\]

This calculation also proves a physical chamber-graph diamond joining the
previous \(D_{55}\) branch. It does not assert global closure.

## Exact local result

All 33 primitive walls are singleton ambient wall orbits. Restriction to the
fixed ten-dimensional \(\mathbb Q(\sqrt5)\)-eigenspace gives

\[
\boxed{
33\text{ ambient walls}
\longrightarrow24\text{ restricted rays}
\longrightarrow12\text{ facets}.}
\]

Exact timelike relative-interior points prove all 12 facets. The remaining
12 rays have exact strictly positive conic expressions in the facet rays.
The facets split into eight root facets and four internal facets.

The root classes are

\[
\begin{gathered}
\{(1,0),(6,0)\},\quad
\{(2,0),(9,0)\},\quad
\{(4,0),(17,0)\},\quad
\{(5,0),(12,0)\},\\
\{(7,0),(16,0)\},\quad
\{(8,0),(10,0)\},\quad
\{(11,0),(15,0)\},\quad
\{(13,0),(18,0)\}.
\end{gathered}
\]

The internal classes are

\[
\boxed{
\{(19,0),(31,0)\},\quad (20,0),\quad (21,0),\quad (23,0).}
\]

The first class is a single eigenspace facet supported by two distinct
ambient walls: the second restricted functional is \(\phi\) times the first.

## Divisor/root interpretation

The ambient chamber has 19 `_minustwo` walls. Exact recovery in the
Néron--Severi lattice and the golden field structure gives:

- 18 integral field-reflection roots;
- nine roots of norm \(2\) and nine of norm \(4-2\phi\);
- field-span rank 9; and
- the unique non-field root `(14,0)`, which is redundant.

All 16 ambient roots supporting the eight root facets are field roots.

## Three simple internal crossings

Exact opposite-side points give

| \(D_{27}\) wall | target representative | target walls | physical identification |
|---:|---:|---:|---|
| `(20,0)` | \(D_{86}\) | 35 | the same marked \(D_{86}\) entered from |
| `(21,0)` | \(D_{14}\) | 40 | the same marked \(D_{14}\) reached from \(D_{55}\) |
| `(23,0)` | \(D_{61}\) | 32 | new on the explored frontier |

The first equality is proved from the cumulative row isometry and complete
transported wall set. The second is the same comparison against the
\(D_{55}\to D_{14}\) state. Thus there is an exact marked-chamber diamond

\[
\require{AMScd}
\begin{CD}
D_{86} @>>> D_{55}\\
@VVV       @VVV\\
D_{27} @>>> D_{14}.
\end{CD}
\]

This is a physical equality, not an inference from the repeated label 14.

## The two-wall simultaneous facet

An infinitesimal move across \(\{(19,0),(31,0)\}\) leaves the initial
Conway chamber through both ambient walls. Crossing either one first gives
one of two exact four-step paths:

\[
\begin{aligned}
&(19,0),(24,0),(20,1),(26,0),
&&27,11,11,74,27,\\
&(31,0),(25,0),(24,0),(16,1),
&&27,74,11,11,27.
\end{aligned}
\]

At every intermediate state the face point is timelike, has exactly two zero
wall pairings, and has no negative pairing. Both paths end at exactly the
same cumulative row isometry and physical wall set. An exact point at scale
\(1/3\) beyond the face is strictly inside that endpoint.

The endpoint again has representative label 27, but it is a different
physical \(D_{27}\): its wall-set fingerprints are

```text
initial  8272b57f633c4a6d3a0474cf53a5c5e279b897d7a1d3d9e513d988bb88935735
opposite 2008055693fd68ed2a0c5ac7778e18edda27c83be5d71c56dbea7b92168b663b
```

This is precisely why representative-index equality is not a valid graph
closure test.

## Reproduction

The external Shimada supplement has SHA-256

```text
fa39e5b59f4a4b2b1b7fe7632e5737d1734ee3eb7e7023ae26bf2b83387afb68
```

Run

```sh
PYTHONPATH=src:. python3 scripts/analyze_shimada_disc11_d27.py /path/to/disc11.txt
SHIMADA_DISC11=/path/to/disc11.txt PYTHONPATH=src:. \
  pytest -q tests/test_shimada_disc11_d27.py
```
