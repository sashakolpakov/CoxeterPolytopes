# Exact completion of the degree-29 four-black stratum

**Status (2026-08-10): every degree-`29` marked passport with at most four
black vertices is impossible; global degree-`29` sharpness remains open.**

Let

\[
 \lambda=\frac{7+2\sqrt5}{29},\qquad
 29\lambda^2-14\lambda+1=0.
\]

As throughout the passport ledger, replace `B` by `1-B` if necessary so
that `B(lambda)=0`; "black" therefore means the fibre containing the marked
point.

The preceding four-black reductions left exactly `3,458` systems in which
one normalized source endpoint lies in the black fibre and the other is a
simple point of the white fibre.  They represent `10,374` marked passports.
The calculation below eliminates all of them.

> **Four-black completion theorem.**  There is no degree-`29` Shabat
> polynomial with the marked black point `lambda`, exactly four distinct
> black roots, and one normalized endpoint in each fibre.

The equations do not impose simplicity on the selected white endpoint, so
the theorem is slightly stronger than the remaining ledger block.  Together
with the independently certified `BB`, white-critical `BW/WB`, and all-status
`WW` no-go results, it closes the complete four-black stratum.  Combining it
with the earlier one-, two-, and three-black no-go theorems gives:

> **Corollary.**  Any degree-`29` Shabat polynomial witnessing the lower bound
> for the `m=5` Gram--Belyi label must have at least five distinct black
> vertices.

This is not a degree-`29` no-go: the strata with at least five black vertices
remain open.

## 1. The `1,729` canonical systems

Write a black-degree record as `(a,b,c,d)`, where `a` is the multiplicity of
the marked root `lambda`, `b` is the multiplicity at the black endpoint, and
the unlabelled remaining degrees satisfy `c<=d`.  There are exactly `1,729`
such records.

For orientation `BW`, scale by `y=x/lambda`.  For orientation `WB`, first
reflect by `x -> 1-x` and then scale by `y=(1-x)/(1-lambda)`.  In either case
the black endpoint is `0`, the marked root is `1`, and

\[
                  G(y)=y^b(y-1)^a(y-u)^c(y-v)^d.
\tag{1}
\]

The other, white endpoint has the fixed coordinate

\[
            z_{BW}=\frac1\lambda,
            \qquad
            z_{WB}=\frac1{1-\lambda}.
\tag{2}
\]

Thus one critical configuration covers both orientations.  The `1,729`
degree records give `2*1,729=3,458` oriented systems.

## 2. Critical-value and endpoint equations

Remove the four forced black-root factors from `G'`, leaving the residual
cubic `D`.  Reduce `G` modulo `D`:

\[
                   G(y)\equiv R_2(u,v)y^2+R_1(u,v)y+R_0(u,v)\pmod D.
\tag{3}
\]

All non-root critical points have one value if and only if

\[
                         R_1(u,v)=R_2(u,v)=0.
\tag{4}
\]

This is scheme-theoretic and includes repeated residual critical points.
For every record the two equations have total degrees `28,27`, bidegrees
`(27,27)`, and

\[
                  \deg_u\operatorname{Res}_v(R_1,R_2)=756.
\tag{5}
\]

At a solution of (4), the common white critical value is `R_0`.  The selected
endpoint belongs to that white fibre only if

\[
                             G(z_{BW})-R_0=0
                 \quad\hbox{or}\quad
                             G(z_{WB})-R_0=0.
\tag{6}
\]

No derivative or simplicity condition at `z` is used.

## 3. Exact extension-field census

For each degree record, choose a split prime for `Q(sqrt(5))` and form

\[
                         A(u)=\operatorname{Res}_v(R_1,R_2).
\tag{7}
\]

Factor `A` over `F_p`.  For every irreducible factor `f(u)`, pass to

\[
                         K_f=\mathbf F_p[u]/(f)
\tag{8}
\]

and compute the full specialized gcd

\[
                         C_f(v)=\gcd_v(R_1,R_2).
\tag{9}
\]

For each irreducible factor `g(v)` of `C_f`, work in the exact tower field

\[
                         L_{f,g}=K_f[v]/(g).
\tag{10}
\]

Both elements `G(z_BW)-R_0` and `G(z_WB)-R_0` are nonzero in every one of the
`28,159` tower orbits.  This checks algebraic-closure points rather than only
`F_p`-rational points.  The first successful prime for each record was:

| split prime | records separating both orientations |
|---:|---:|
| 101 | 1,529 |
| 109 | 174 |
| 131 | 23 |
| 139 | 2 |
| 149 | 1 |

## 4. Why reduction cannot lose a characteristic-zero solution

The characteristic-zero degree bound in (5) is `756`, and every selected
finite-field resultant retains degree exactly `756`.  After division by its
unit leading coefficient, it is monic, so a putative `u` is integral at the
selected prime.  On every resultant component, at least one specialized
critical equation retains degree `27` in `v`; its unit leading coefficient
makes a common `v` integral.  Finally, both coordinates in (2) are local
units because the selected primes avoid `2` and `29`.

Consequently, every characteristic-zero solution of (4) and (6) would reduce
to one of the exact tower components checked above.  Since the relevant
endpoint element is nonzero on every component, no such solution exists.

## 5. Exact boundary after completion

The final block contributes three simple-white marked passports for each
oriented degree system, one from each possible white partition.  Hence the
`3,458` systems account for `3*3,458=10,374` marked records.  The complete
boundary is

```text
four_black_status = IMPOSSIBLE
four_black_systems_remaining = 0
four_black_marked_passports_remaining = 0
degree_29_black_vertex_lower_bound = 5
global_marked_passports_eliminated = 48,343
global_marked_passports_remaining = 42,562,335
global_degree_29_status = OPEN
degree_29_witness_found = false
complete_degree_29_no_go_claimed = false
```

## 6. Reproducibility

The exact producer and canonical archive are:

- [`certify_gram_belyi_m5_degree29_four_black_one_black_simple.py`](../scripts/certify_gram_belyi_m5_degree29_four_black_one_black_simple.py);
- [`gram_belyi_m5_degree29_four_black_one_black_simple.json`](../data/enumerations/gram_belyi_m5_degree29_four_black_one_black_simple.json).

The archive stores every monic degree-`756` resultant, every specialized
component degree, and both endpoint nonvanishing results for all tower
orbits.  The replay is a bounded `1,729`-case exact census and launches no
plane-tree search.  No timeout or failed search is used as evidence.
