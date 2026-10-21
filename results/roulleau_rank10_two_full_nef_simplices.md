# Two rank-ten K3 full-nef Coxeter simplices: `T9` and `Q9`

**Status date:** 2026-08-09.  **Verdict:** two positive ledger rows,
finite-volume and noncompact, with exact full-nef K3 realizations.  They are
known abstract simplices and published K3 configurations; the new item is the
explicit project-level same-polytope alignment and literal certificate, not a
worldwide-priority claim.

## Result

Two rows of Roulleau's nine-member Picard-rank-ten list had not been expanded
as individual positive entries in the project ledger:

| ID | exact Picard lattice | complete walls | exact simplex | finite / ideal vertices | status |
|---|---|---:|---|---:|---|
| `R10-UE8-NEF-H9-10` | `U+E8(-1)` | 10 irreducible `(-2)`-curves | JKRT `T9` | `9 / 1` | finite-volume, noncompact |
| `R10-UD8-NEF-H9-10` | `U+D8(-1)` | 10 irreducible `(-2)`-curves | JKRT `Q9` | `7 / 3` | finite-volume, noncompact |

For either row, let `A1,...,A10` be Roulleau's curves and put

\[
 G_{ii}=-2,\qquad G_{ij}=1\text{ on the listed edges},\qquad
 G_{ij}=0\text{ otherwise}.
\]

The edge sets, in one-based notation, are

```text
U+E8:  12 23 34 45 56 67 78 79 9-10
U+D8:  13 23 34 45 56 67 78 79 9-10.
```

Thus the first diagram is literally `T9`.  Swapping labels 2 and 3 in the
second makes it literally `Q9`.  Both diagrams have only `m=2,3` relations;
in particular, `Q9` must not be conflated with the different `S9/P9` simplex,
which has an `m=4` edge.

The exact matrix invariants are

```text
             det(G)        Smith(G)                 inertia
U+E8           -1          1^10                     (1,9,0)
U+D8           -4          1^8,2,2                  (1,9,0).
```

The full matrices and their canonical hashes are frozen in the
[certificate](../data/enumerations/roulleau_rank10_two_full_nef_simplices.json).

## The algebraic surfaces and exact cone equality

Roulleau proves that the ten displayed curves generate the stated
Néron--Severi lattice and are the complete irreducible `(-2)`-curve
configuration.  He also supplies ample divisors in the wall basis:

```text
U+E8: (30,61,93,126,160,195,231,115,153,76), square 1240;
U+D8: (15,15,31,33,36,40,45,22,29,14),        square  280.
```

Exact multiplication gives `G D=(1,...,1)` in both cases.  Hence each
published divisor lies strictly inside the chosen chamber and fixes the time
sheet.  The K3 ample-cone theorem now gives the strong conclusion

\[
 \operatorname{Nef}(X)=
 \{x\in\overline{\mathcal C}^{+}:x\cdot A_i\geq0,
   \ 1\leq i\leq10\}.
\]

This is equality with the **entire** nef cone, not a divisor subsystem or an
automorphism-fixed slice.  No surface automorphism is used.

The surfaces themselves have concrete algebraic models.  In the `E8` row,
the elliptic fibre is

\[
 A_2+2A_3+3A_4+4A_5+5A_6+6A_7+3A_8+4A_9+2A_{10},
\]

with section `A1`; `|4F+2A1|` realizes the K3 as a double cover of
`F4`, branched over its negative section and a member of `|3s+12f|` with an
`E8` singularity.  In the `D8` row the fibre is

\[
 A_1+A_2+2(A_3+\cdots+A_7)+A_8+A_9,
\]

with section `A10`; the analogous branch curve has a `D8` singularity.
These models and complete curve graphs are in Roulleau's
[*Atlas of K3 surfaces with finite automorphism group*](https://arxiv.org/abs/2003.08985),
rank-ten section.

## Independent cone and Vinberg gates

Since each wall Gram is nonsingular, the ten-wall chamber in `H9` is a
simplex.  The columns of `G^{-1}` are its opposite projective vertex rays.
Their exact squared norms are

```text
T9:  0,2,6,12,20,30,42,10,18,4;
Q9:  0,0,2,3,4,5,6,1,2,0.
```

There are therefore no spacelike/hyperideal rays: `T9` has nine finite
vertices and one cusp, while `Q9` has seven finite vertices and three cusps.
For every wall `i`, the sum of the other nine dual rays has square greater
than zero, pairs `0` with wall `i`, and pairs `1` with every other wall.  This
gives ten exact strict relative-facet witnesses and proves irredundancy.

The generic literal Vinberg verifier separately enumerates all 45 elliptic
rank-eight edges.  Its output is

| simplex | edge endpoint types `(finite,ideal)` | connected affine components | bad edges | bad full parabolic neighbourhoods |
|---|---|---:|---:|---:|
| `T9` | `36*(2,0)+9*(1,1)` | 1 | 0 | 0 |
| `Q9` | `21*(2,0)+21*(1,1)+3*(0,2)` | 3 | 0 | 0 |

Consequently both chambers satisfy Vinberg's complete edge-endpoint and
parabolic-neighbourhood tests.  They are finite-volume **and** noncompact;
the cusps are not being mistaken for evidence of infinite volume.

## Lineage and denominator

The abstract simplex diagrams are classical.  The exact graph names `T9`
and `Q9` occur in Johnson--Kellerhals--Ratcliffe--Tschantz,
[*The size of a hyperbolic Coxeter simplex*](https://doi.org/10.1007/BF01238563),
which computes the covolumes of the classified hyperbolic Coxeter simplices.
The two CoxIter graph snapshots are hash-pinned in the certificate and agree
with the K3 graphs entry by entry under the permutations above.

The comparison with the `118`-type theorem is deliberately scoped.  Before
this batch, the global theorem covered these rows collectively, but the
dedicated positive packages expanded the eight compact rows and the nineteen
rank-11--19 rows.  Roulleau has nine rank-ten types; this batch expands two.
The exact finite batch of the other seven is now also complete:

```text
U+E7+A1, U+D4^2, U+D6+2A1, U(2)+D4^2,
U+D4+4A1, U+8A1, U+A2+E6.
```

All seven have exact full-nef and literal finite-volume certificates in the
[rank-ten completion theorem](roulleau_rank10_remaining_full_nef.md).  This is
a row-level denominator closure, not a claim that no unrelated construction
elsewhere in the large project ledger is isometric to one of these chambers.

## Replay

```bash
PYTHONPATH=.:src python3 scripts/certify_roulleau_rank10_two_full_nef_simplices.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_roulleau_rank10_two_full_nef_simplices.py
```

The first command recompiles the exact C++ Vinberg verifier and rechecks the
entire archived payload.
