# Independent hostile audit: silver prime-seven `D24+` through `H = 500`

## Verdict

**PASS.**  A clean-room exact replay found no defect in the frozen
`H <= 500` census or its infinite-volume certificate for

\[
  \langle -(2\sqrt2-1)\rangle\oplus D_{24}^{+}.
\]

The audit independently reproduces all 64 admissible controls through height
500, the 16 controls above the pinned `H = 300` frontier, all new square
shells, all 4,026 newly accepted walls, and the final count of 4,374 walls.
It also independently confirms that the displayed replacement ray is a
rank-24 hyperideal extreme ray.  The conclusion remains deliberately local:
the complete `H <= 500` prefix has infinite volume, while the eventual
Vinberg chamber is unresolved.

## Clean-room boundary

`scripts/independent_hostile_audit_silver_p7_d24_plus_h500.py` imports neither
the production search nor the project's number-field implementation.  It
uses standalone `Fraction` pairs for `Q(sqrt(2))`, constructs an explicit
integral `D24+` basis, verifies its even unimodular Gram matrix, and enumerates
the root residue and all three nonzero spinor residues directly in
orthonormal coordinates.

The local norm argument was checked separately.  With

```text
pi = sqrt(2),  epsilon = 1+sqrt(2),  delta = 3+sqrt(2),
a = 2sqrt(2)-1,  epsilon*a = delta,  N(delta)=7,
```

reflection integrality excludes prime support away from `(pi)` and
`(delta)`.  At `(delta)`, unimodularity forces the spatial part into
`delta*D24+`, and primitivity makes the norm valuation exactly one.  At
`(pi)`, a unit control gives valuation zero; otherwise evenness,
unimodularity, and primitivity force valuation exactly two.  Since
`Z[sqrt(2)]` has class number one and each totally positive unit is a square,
the four canonical norms are exactly

```text
1, 2, 3+sqrt(2), 6+2sqrt(2).
```

## Completeness replay

For a control `x=m+n*sqrt(2)`, the height and conjugate-target inequalities
give an exact ellipse

```text
m^2 + 2n^2 <= ceil((x_0^2 bound + x_1^2 bound)/2).
```

The four ceilings, in canonical norm order, are `137, 274, 604, 1208`.
Exhausting those boxes gives exactly 64 controls through `H = 500`, of which
16 lie in `(300,500]`; their complete exact set equals the frozen set.

Spatial enumeration uses the analogous trace bound.  For an integral
coordinate `z=m+n*sqrt(2)` in a target `T=A+B*sqrt(2)`, completeness follows
from `m^2+2n^2 <= A`.  For a half-spinor coordinate the doubled coefficients
satisfy `m^2+2n^2 <= 4A`.  Exact two-embedding remainder tests, parity tests,
and a 24-slot multiplicity recursion then exhaust every dominant `D24+`
square partition.

At the dense height `181+128sqrt(2)`, the two control levels give:

| Norm | Raw reflective | Imprimitive `sqrt(2)` | Primitive |
|---|---:|---:|---:|
| `1` | 83 | 0 | 83 |
| `2` | 4,090 | 83 | 4,007 |
| **packet** | **4,173** | **83** | **4,090** |

Independent filtering against the 348 pinned earlier walls retains exactly
4,003 roots, with exact set equality against the frozen archive.  All
8,010,003 unordered pairs among those 4,003 survivors are nonpositive.
Consequently they are a genuinely simultaneous packet; insertion order did
not suppress or manufacture a wall.  Including the later 6- and 17-wall
packets gives 8,010,154 unordered same-height checks in total.

The integer matrix replay evaluates 17,587,030 exact entries.  Its independent
worst-case intermediate bound is 147,376, versus `2^63-1`; the largest actual
cleared pairing coefficient is 4,632.  Standalone scalar arithmetic replays
76,768 filtering boundaries, another 1,490 same-height zero pairs, three
strict inequalities for every retained root, and one positive witness for
each rejected candidate.

## Rays and exact geometry

The old `H = 300` replacement ray is cut by the extension exactly as reported:
3,952 negative, zero zero, and 74 positive pairings.  Its largest new pairing
is

\[
  \frac{594-411\sqrt2}{7}>0.
\]

For

\[
  c=4-\frac{12}{5}\sqrt2,\qquad
  v=(1;c,c,c,c,c,0^{19}),
\]

clean-room scalar pairing against all 4,374 walls gives 4,350 negative, 24
zero, and no positive values.  The zero facets are

```text
0,1,2,3,5,...,23,4361.
```

Exact Gaussian elimination over `Q(sqrt(2))` gives rank 24 for their wall
functionals.  Their common nullspace is therefore the ray spanned by `v`, so
feasibility makes it an extreme ray.  Its norm is

\[
  (v,v)=\frac{693}{5}-98\sqrt2>0,
\]

making it hyperideal.  Facet 4,361 independently remains the limiting wall
for the five-equal-coordinate construction.

## Reproducibility

```text
python3 scripts/independent_hostile_audit_silver_p7_d24_plus_h500.py
pytest -q tests/test_silver_p7_d24_plus_h500_independent_hostile_audit.py
```

The machine-readable certificate is
`data/enumerations/silver_p7_d24_plus_h500_independent_hostile_audit.json`.
The independent replay intentionally trusts only the already pinned
`H = 300` wall stream as its starting frontier; it reconstructs every control,
candidate, acceptance decision, and geometric claim in the `(300,500]`
extension.
