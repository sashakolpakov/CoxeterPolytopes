# Fresh golden H9 lattice modifications

This follow-up changes the lattice rather than extending the earlier
215-wall height-1000 search. It has not produced a compact example. It does
identify two exact duplicate routes before computation and tests four new
quadratic spaces with independent escape certificates.

## The tempting unimodular glue is a duplicate

Let `O=Z[phi]` and `L=<-2phi>+E8+<2>`, with binary basis `(p,u)`.
The vector

`v=(p+phi² u)/2`, with `v²=phi²`,

defines an integral overlattice `M=L+Ov`. In `M`, the basis

`w=(-phi p-u)/2`, `f=((phi-1)p+phi u)/2`

has binary Gram matrix `diag(-phi,1)` and coordinate determinant `-1/2`.
Thus `M` is unimodular and isometric to `<-phi>+E8+<1>`.

This does **not** yield a new reflection group. Over `O/(2)=F4`, an integral
overlattice of the binary block requires `b²=phi a²`; this has one
projective solution, so `M` is the unique proper integral overlattice.
Every isometry of `L` therefore preserves `M`. Conversely,
`L={x in M : x² in 2O}`, so every isometry of `M` preserves `L`.
Consequently the two integral orthogonal groups, including their full
reflection subgroups, are identical in the common quadratic space. In
particular this glue cannot rescue the previously unclosed full chamber.

A second tempting variant is also a duplicate. In the binary form
`diag(-phi,2)`, the vectors

`w=(2phi,phi²)`, `f=(phi,phi)`

have Gram `diag(-2phi,1)` and determinant one. Hence
`<-phi>+E8+<2>` is already the previously searched odd-A1 host.

## New determinant hosts

The fresh hosts are `L_alpha=<-alpha>+E8+<1>`, for

`alpha=sqrt5=2phi-1`, `3phi`, `1+4phi`, and `-1+3phi`.

All have signature `(9,1)` at the chosen real place and positive-definite
conjugate. Their determinant ideals have odd valuation at, respectively,
the prime over5, the inert prime3, and the two different primes over11.
These distinguish their determinant square classes from one another and
from the earlier natural and odd-A1 hosts. In rank ten the determinant
square class is also invariant under scalar similarity, so these are
genuinely different quadratic spaces.

For a primitive integral reflection, its norm `q` divides `2 det(B)`.
The search lists every canonical totally positive divisor of `2 det(B)`;
unit-square normalization is valid because every totally positive unit of
`O` is a square. The resulting complete **candidate** alphabets are:

| alpha | norm candidates |
|---|---|
| `2phi-1` | `1,2,2+phi,4+2phi` |
| `3phi` | `1,2,3,6` |
| `1+4phi` | `1,2,3+phi,6+2phi` |
| `-1+3phi` | `1,2,4-phi,8-2phi` |

Empty classes are harmlessly retained. For example, modulo the inert
prime3, the norm-three and norm-six equations in the `3phi` host force a
nonsquare in `F9` to be a square; those classes have no primitive mirrors.

The height-zero arrangement is fully known for every host: `E8+A1`, with
the A1 root having norm one. Its spatial lattice is unimodular, so a
primitive spatial reflection has normalized norm one or two. Writing an
E8 coordinate as `a+phi b` gives trace norm
`2||a+b/2||²+(5/2)||b||²`. Its minimum is four, with equality only for a
rational E8 root. The trace budgets then force exactly the E8 roots and
the last-coordinate norm-one roots; `u²=2` has no solution in the field.

## Completed search results

The height is `H=alpha*x²/q`. Both-embedding bounds enumerate every possible
control coordinate. Reflection integrality is enforced on the time
coordinate and every E8/last-coordinate pairing. Exact Schur-complement
recursion and chronological nonpositive-pairing gates then produce the
listed prefixes.

| alpha | completed H | walls | finite vertices | hyperideal rays |
|---|---:|---:|---:|---:|
| `2phi-1` | 30 | 16 | 31 | 117 |
| `3phi` | 30 | 13 | 3 | 77 |
| `1+4phi` | 30 | 16 | 3 | 170 |
| `-1+3phi` | 30 | 14 | 24 | 54 |
| `2phi-1` | 100 | 43 | not enumerated | exact escape |

All fully enumerated rows have zero cusps and zero opposite-sheet rays.
The height-100 `sqrt5` continuation took about ten seconds; its exact
spacelike escape has norm `(38486-23762phi)/841 > 0`.
The fourth determinant variant is recorded in
`fresh_golden_neighbor_m1p3phi_h30.json`.

The standalone verifier `scripts/fresh_golden_neighbor_verify.py` imports
no project arithmetic or search code. Python `Fraction` arithmetic checks
every root norm, every reflection coefficient, all pairwise acuteness,
height ordering, the strict timelike interior, and every inequality of the
spacelike witness. The four first-height certificates and the height-100
certificate pass.

These certificates prove that the displayed finite prefixes have infinite
volume. They are not global nonreflectivity proofs and do not exclude later
compact completion.

## Reproduction

```
PYTHONPATH=src:. python3 scripts/fresh_golden_neighbor_search.py \
  --alpha=-1,2 --height 100 \
  --output results/fresh_golden_neighbor_sqrt5_h100.json
PYTHONPATH=src:. python3 scripts/fresh_golden_neighbor_geometry.py \
  results/fresh_golden_neighbor_sqrt5_h100.json --escape-only
python3 scripts/fresh_golden_neighbor_verify.py \
  results/fresh_golden_neighbor_sqrt5_h100.json
PYTHONPATH=src:. python3 -m pytest -q tests/test_fresh_golden_neighbor.py
```

The generic search exactly reproduces the earlier odd-A1 ordered 17-wall
prefix, including its primitive normal scaling. Eight tests pass, covering
that calibration, both explicit lattice identifications, the unique F4
gluing line, and five independent certificate replays. Search checkpoints
are written after every completed height packet, so interruption preserves
an accurately labeled bounded result.
