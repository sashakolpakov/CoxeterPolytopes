# Independent hostile audit: silver `DualExtremal(28,3)`

**PASS.** This replay never imports the production owner module. It parses the
Gram literal as inert source data and independently reconstructs every claim
needed for the bounded prefix through `H=93+66sqrt(2)`.

The 28-dimensional Gram is even integral symmetric, has determinant `3^26`,
Smith form `1^2,3^26`, rank two modulo three, minimum 8, and 1,260 signed
minimal vectors. An independent GF(3)-kernel lift and column-HNF computation
gives the pairing-radical preimage index 9; PARI finds its exact minimum 12
and 65,520 signed minimal vectors.

The local replay gives exactly eight dyadic/prime-seven/discriminant-three
norm classes. It enumerates 51 controls on 35 heights. Every earlier control
has an empty rigorous necessary Gram superset. The endpoint contains exactly
six types:

| norm | time | Gram `(u²,v²,u.v)` | augmented maximum | `d/2` margin |
|---:|---:|---:|---:|---:|
| 3 | 9 + 6*sqrt(2) | (138, 72, 99) | 24 + 18*sqrt(2) | 3 + 3*sqrt(2) |
| 6 | 12 + 9*sqrt(2) | (300, 132, 198) | 45 + 30*sqrt(2) | 3/2 |
| 6 | 12 + 9*sqrt(2) | (288, 138, 198) | 44 + 31*sqrt(2) | 2 - 1/2*sqrt(2) |
| 6 | 12 + 9*sqrt(2) | (276, 144, 198) | 44 + 32*sqrt(2) | 2 - 1*sqrt(2) |
| 6 | 12 + 9*sqrt(2) | (264, 150, 198) | 40 + 31*sqrt(2) | 4 - 1/2*sqrt(2) |
| 6 | 12 + 9*sqrt(2) | (252, 156, 198) | 33 + 27*sqrt(2) | 15/2 + 3/2*sqrt(2) |

Here the witness is explicit: `d=e0`, so `d²=8`, and `z=(1,d/2)` has
`z²=3-2sqrt(2)>0`. All six margins are strictly positive. Therefore the
complete necessary prefix has an open hyperideal patch and infinite volume.

No production defect was found. This remains a bounded negative result: it
does not determine later walls, Vinberg termination, or the eventual chamber.
