# Exact silver prime-seven `D24+` prefix through `H = 500`

## Result

The exact Vinberg prefix for

\[
  \langle -(2\sqrt2-1)\rangle\oplus D_{24}^{+}
\]

was extended from the independently frozen `H = 300` wall list through
algebraic height `H = 500`.  The complete bounded prefix has **4,374 walls**:
24 height-zero `D24` walls and 4,350 accepted positive-height roots.

The `H = 500` prefix still has infinite volume.  An explicit exact
rank-24 hyperideal extreme ray is

\[
 v=(1;c,c,c,c,c,0^{19}),\qquad
 c=4-\frac{12}{5}\sqrt2,
\]

with

\[
 (v,v)=\frac{693}{5}-98\sqrt2>0.
\]

All 4,374 wall pairings with `v` are nonpositive: 4,350 are negative and
24 are zero.  The active wall functionals have exact rank 24.  The unique
positive-height active wall is facet 4,361, in the final nonempty packet at
`H = 248 + 176*sqrt(2)`; the other 23 active facets are height-zero `D24`
walls.  This is a direct exact non-finite-volume certificate and does not
require a full double-description expansion through the 4,003-wall packet.

## Complete bounded census

The local mirror-norm proof remains the frozen four-class proof for norms
`1`, `2`, `3 + sqrt(2)`, and `2*(3 + sqrt(2))`.  Above `H = 300`, all 16 new
admissible controls and every `D24` root/spinor glue class were enumerated.
There are ten new simultaneous-height packets.  Only three are nonempty:

| Height | Raw reflective | Imprimitive `sqrt(2)` | Primitive candidates | Accepted walls | Accepted spinor |
|---|---:|---:|---:|---:|---:|
| `181 + 128*sqrt(2)` | 4,173 | 83 | 4,090 | 4,003 | 1,853 |
| `206 + 146*sqrt(2)` | 7 | 0 | 7 | 6 | 1 |
| `248 + 176*sqrt(2)` | 24 | 0 | 24 | 17 | 3 |

Thus the extension contributes 4,026 walls, taking the accepted-root count
from 324 to 4,350 and the final wall count from 348 to 4,374.

The old `H = 300` replacement ray does not survive: against the 4,026 new
walls it has 3,952 negative, zero zero, and 74 positive pairings.  Its largest
new pairing is

\[
 \frac{594}{7}-\frac{411}{7}\sqrt2>0.
\]

## Exact dense-packet audit

The 4,090 candidates at `H = 181 + 128*sqrt(2)` are filtered and the 4,003
survivors are checked simultaneously, before any member is inserted as a
wall.  Pairings are represented as integer coefficient pairs `A+B*sqrt(2)`
over positive common coordinate denominators.  Matrix products use signed
64-bit integers only after a symbolic absolute-value bound is proved for the
entire Lorentz pairing expression.  Signs are then decided exactly by the
signs of `A`, `B`, and `A^2-2*B^2`; no floating-point approximation enters
selection.

The audit checks 8,010,154 unordered equal-height survivor pairs.  Across
filtering and simultaneous-packet verification it evaluates 17,587,030 exact
matrix entries.  The largest proved intermediate absolute bound is 142,352,
far below `2^63-1`; the largest actual pairing coefficient is 4,632, far below
the explicit exact-sign squaring limit of 2,000,000,000.  The independent
scalar `Silver` arithmetic replays all 76,768 zero boundaries plus 12,078
strict retained inequalities.

The accepted-root stream, per-packet streams, every accepted root vector,
the overflow audit, the old-ray cut audit, and the replacement-ray rank and
pairings are frozen in
`data/enumerations/silver_p7_d24_plus_h500.json`.

## Scope

This result is exhaustive only through algebraic height 500.  It proves that
the complete bounded prefix has infinite volume, not that the eventual
Vinberg chamber has infinite volume: a wall above height 500 may cut the
displayed hyperideal ray.  The eventual chamber remains unresolved.  No
global nonreflectivity or compactness claim is made.
