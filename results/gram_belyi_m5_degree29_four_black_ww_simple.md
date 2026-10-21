# Exact no-go for every four-black `WW` endpoint status

**Status (2026-08-10): all `640` black-degree master systems, subsuming
`2,560` endpoint-status systems and `8,960` marked passports, are impossible;
global degree-`29` sharpness remains open.**

Let

\[
 \lambda=\frac{7+2\sqrt5}{29},\qquad
 29\lambda^2-14\lambda+1=0.
\]

The first four-black reduction identified `640` master systems in which both
normalized endpoints are simple white points.  The equations used below do
not impose simplicity and therefore close the entire `WW` locus:

> **Four-black `WW` theorem.**  No degree-`29` Shabat polynomial with exactly
> four black vertices, one of them marked by `lambda`, can have both
> normalized source points in the white fibre.

Indeed, the calculation allows either selected white point to be critical.
One black-degree master system therefore eliminates all four endpoint-status
systems `SS,SC,CS,CC`.  Across the three white partitions these four statuses
represent `3+4+4+3=14` marked passports per black-degree record.  Hence the
`640` exact certificates eliminate `2,560` endpoint-status systems and
`640*14=8,960` marked records.

## 1. Canonical critical equations

For a marked black-degree record `(a,b,c,d)`, translate the marked black root
to zero and scale one other root to one.  Every configuration is represented
by

\[
                         G(y)=y^a(y-1)^b(y-u)^c(y-v)^d,
 \qquad a+b+c+d=29.
\tag{1}
\]

There are exactly `640` such marked multiplicity records.  If the original
endpoints become `y_0,y_1`, their ratio is fixed:

\[
 q=\frac{y_1}{y_0}=-\frac{1-\lambda}{\lambda}
   =29\lambda-13,\qquad q^2+12q+16=0.
\tag{2}
\]

Remove the forced black-root factors from `G'`, obtaining the residual cubic
`D`.  Write

\[
                         G(y)\equiv R_2(u,v)y^2+R_1(u,v)y+R_0(u,v)pmod D.
\tag{3}
\]

All non-root critical points have one value if and only if

\[
                         R_1(u,v)=R_2(u,v)=0.
\tag{4}
\]

This includes every repeated-root factorization of `D`.  Homogeneity gives

\[
 \deg R_1=28,\qquad \deg R_2=27.
\tag{5}
\]

Both equations have degree at most `27` in either one of `u,v`.  For `R_1`,
the apparent pure degree-`28` term vanishes: coalescing the other three black
roots at zero makes both `G` and `D` divisible by `y^2`, so the remainder has
only a `y^2` term.  Every selected reduction has exact bidegrees `(27,27)`.

For a critical configuration put

\[
                         H_{u,v}(y)=G(y)-R_0(u,v).
\tag{6}
\]

The two normalized endpoints could lie in the common white fibre only if

\[
                         H_{u,v}(y)=H_{u,v}(qy)=0
\tag{7}
\]

for some `y`.

## 2. Exact resultant and extension-field census

At a split prime form

\[
                         A(u)=\operatorname{Res}_v(R_1,R_2).
\tag{8}
\]

Bezout gives `deg A<=28*27=756`; every certificate has exact degree `756`.
Factor `A` over `F_p`.  For each irreducible factor `f(u)`, work in the exact
field

\[
                         K_f=\mathbf F_p[u]/(f)
\tag{9}
\]

and compute the full polynomial gcd

\[
                         C_f(v)=\gcd_v(R_1,R_2).
\tag{10}
\]

Factors caused only by a resultant collision at infinity have `C_f=1` and
are retained in the ledger rather than silently discarded.  For every
nonconstant irreducible factor `g(v)` of `C_f`, the calculation then works in
the tower field

\[
                         L_{f,g}=K_f[v]/(g)
\tag{11}
\]

and computes

\[
                  \gcd_y\bigl(H_{u,v}(y),H_{u,v}(qy)\bigr).
\tag{12}
\]

Every one of the `10,696` irreducible `(u,v)` orbits has gcd degree zero.
This is an algebraic-closure test in exact finite extensions, not a search for
`F_p`-rational points.

The witness-prime histogram is

| split prime | systems separated |
|---:|---:|
| 101 | 572 |
| 109 | 58 |
| 131 | 10 |

## 3. Why reduction cannot lose a solution at infinity

Three monicity checks make the finite-field result a characteristic-zero
certificate.

1. The characteristic-zero degree of `A` is at most `756`, while the selected
   reduction retains degree `756`.  Its leading coefficient is therefore a
   local unit; after scaling `A` is monic, so any root `u` is integral.
2. On every relevant factor of `A`, at least one of `R_1,R_2` retains degree
   `27` in `v`.  Its leading coefficient is a unit in the corresponding
   residue field, so a common root `v` is integral.
3. `H_{u,v}(y)` is monic of degree `29`, so a common endpoint coordinate `y`
   is integral.

Thus any characteristic-zero solution would reduce to one of the exact tower
components checked in (12).  Since every tower gcd is one, no such solution
exists.

## 4. Updated exact boundary

The four-black ledger now has only one block:

| remaining block | systems | marked passports |
|---|---:|---:|
| one black endpoint, other white simple (`BW+WB`) | 3,458 | 10,374 |
| **total** | **3,458** | **10,374** |

Together with all preceding reductions, `37,969` of the global `42,610,678`
marked degree-`29` records are excluded.  The exact global remainder is
`42,572,709`.  Therefore

```text
four_black_WW_status = IMPOSSIBLE
four_black_status = PARTIAL
global_degree_29_status = OPEN
degree_29_witness_found = false
complete_four_black_no_go_claimed = false
complete_degree_29_no_go_claimed = false
```

## 5. Reproducibility

The exact producer and archive are:

- [`certify_gram_belyi_m5_degree29_four_black_ww_simple.py`](../scripts/certify_gram_belyi_m5_degree29_four_black_ww_simple.py);
- [`gram_belyi_m5_degree29_four_black_ww_simple.json`](../data/enumerations/gram_belyi_m5_degree29_four_black_ww_simple.json).

The archive retains the monic degree-`756` critical resultants and every
component/orbit degree.  A replay is casewise and bounded, launches no
plane-tree census, and uses no timeout as evidence.  The only remaining
four-black target is the `3,458`-system `BW+WB` block with a simple white
endpoint (`1,729` systems per orientation).
