# Exact quantified reduction of the degree-29 four-black stratum

**Status (2026-08-10): `23,660` of `42,994` marked four-black records are
eliminated exactly; degree-`29` sharpness remains open.**

Let

\[
 \lambda=\frac{7+2\sqrt5}{29},\qquad
 29\lambda^2-14\lambda+1=0.
\]

The preceding three-black theorem implies that any degree-`29` witness,
after complementing so that `B(lambda)=0`, has at least four points in that
black fibre.  The present bounded calculation attacks the entire four-black
stratum at the level of coefficient systems, without enumerating plane trees.
It proves:

> **Four-black partial theorem.**  A degree-`29` witness with exactly four
> black vertices cannot have both normalized endpoints black.  Nor can it
> have exactly one black endpoint while the other endpoint is a non-root
> critical point.  Therefore a surviving four-black witness must have either
> one black and one *simple* white endpoint, or two white endpoints.

This eliminates `6,734` of the `12,752` distinct necessary coefficient
systems, representing `23,660` marked passport records.  Exactly `6,018`
systems representing `19,334` four-black records remain.

## 1. Exact passport/system quotient

The white partition has length `26`, so its excess ramification is three.
There are exactly three possibilities:

```text
(4,1^25),       (3,2,1^24),       (2,2,2,1^23).
```

For one marked white endpoint these partitions have respectively `2,3,2`
possible ramification degrees.  For two ordered white endpoints they have
`3,7,4` choices.  On the black side there are:

- `640` multiplicity records with the `lambda`-part distinguished;
- `1,729` records after additionally choosing one black endpoint degree;
- `3,276` records after choosing two ordered black endpoint degrees.

These reproduce the raw marked-passport census:

| endpoint colors | factorization | marked records |
|---|---:|---:|
| `BB` | `3276 * 3` | 9,828 |
| `BW` | `1729 * 7` | 12,103 |
| `WB` | `1729 * 7` | 12,103 |
| `WW` | `640 * 14` | 8,960 |
| total |  | 42,994 |

For coefficient elimination, the exact white ramification index can be
forgotten while retaining whether an endpoint is simple (`S`) or critical
(`C`).  The common-critical-value equations themselves determine the
factorization type of the residual derivative cubic.  The endpoint-status
counts across the three white partitions are

| status | marked choices across partitions |
|---|---:|
| one endpoint `S` / `C` | 3 / 4 |
| ordered `SS,SC,CS,CC` | 3,4,4,3 |

Thus the `42,994` marked passports map to only

\[
 3276+2(2\cdot1729)+4\cdot640=12,752
\tag{1}
\]

necessary master systems.  This quotient is used only for no-go proofs.  It
is not a claim that every master system constructs every exact multiplicity
passport represented by it; exact ramification indices would add
multiplicity and inequation conditions on the residual cubic.

## 2. All `BB` systems are impossible

For two black endpoints write

\[
 F(x)=(x-\lambda)^a x^b(x-1)^c(x-t)^d,
 \qquad a+b+c+d=29.
\tag{2}
\]

There are exactly `binom(28,3)=3,276` positive ordered exponent quadruples,
one for every coefficient system.  Remove the forced black-root factors from
`F'`; the residual factor `D` is cubic.  Write

\[
                         F(x)\equiv R_2(t)x^2+R_1(t)x+R_0(t)\pmod D.
\tag{3}
\]

All non-root critical points have one value if and only if

\[
                         R_2(t)=R_1(t)=0.
\tag{4}
\]

This remains true when `D` has multiple roots.  Under simultaneous scaling
of all black roots, `R_j` is homogeneous of degree `29-j`.  Hence
`deg_t R_2<=27`.  The apparent degree-`28` leading term of `R_1` is zero:
setting the other three roots to zero makes both `F` and `D` divisible by
`x^2`, and the remainder has only an `x^2` term.  Thus also
`deg_t R_1<=27`.

Every modular certificate retains degree `27` in both equations and has
full polynomial gcd one.  The split-prime histogram is

| prime | systems separated |
|---:|---:|
| 101 | 3,156 |
| 109 | 106 |
| 131 | 14 |

Degree preservation makes these exact nonzero-resultant certificates over
`Q(lambda)`, not finite-field rational-point searches.  Since one master
system covers each of the three white partitions, this eliminates all
`9,828` marked `BB` records.

## 3. One black and one critical white endpoint

Reflection `x -> 1-x` lets both endpoint orientations be treated in the form

\[
 F(x)=(x-\alpha)^a x^b(x-r)^c(x-s)^d,
 \qquad \alpha\in\{\lambda,1-\lambda\},
\tag{5}
\]

where the white endpoint `x=1` is critical.  If `D` is the residual cubic,
then

\[
                         D(1)=C(r)s+K(r),
\tag{6}
\]

with `C,K` linear.  Their linear resultant is nonzero at every selected
certificate prime, excluding the only branch lost by setting `s=-K/C`.

Let `R_1,R_2` again denote the two nonconstant coefficients of `F mod D`.
Homogeneity gives the cleared univariate equations

\[
 \begin{aligned}
 P_1(r)&=C(r)^{28}R_1\bigl(r,-K(r)/C(r)\bigr),\\
 P_2(r)&=C(r)^{27}R_2\bigl(r,-K(r)/C(r)\bigr).
 \end{aligned}
\tag{7}
\]

Their degree bounds are `55` and `54`.  For `P_1`, the only possible
degree-`56` contribution would come from the pure `r^28` coefficient of
`R_1`; the same three-root coalescence argument used above makes it zero.

At each split prime, `56` ordinary values determine `P_1,P_2` exactly by
Newton interpolation, and the certificate verifies the identities at every
remaining residue value with `C(r)!=0`.  It then requires exact degrees
`55,54` and full gcd one.  This excludes all `1,729` systems in each endpoint
orientation.  Their combined prime histogram is

| prime | systems separated |
|---:|---:|
| 101 | 3,285 |
| 109 | 164 |
| 131 | 9 |

A critical white endpoint has `4` marked ramification choices across the
three white partitions.  The `3,458` eliminated systems therefore represent

\[
                         2\cdot1729\cdot4=13,832
\tag{8}
\]

marked passports.

## 4. Exact remaining boundary

The four-black remainder is:

| block | systems | marked passports |
|---|---:|---:|
| one black, other white simple (`BW+WB`) | 3,458 | 10,374 |
| two white, at least one critical | 1,920 | 7,040 |
| two white, both simple | 640 | 1,920 |
| **total** | **6,018** | **19,334** |

Combining this result with the completed one-, two-, and three-black no-goes
eliminates `29,009` of the global `42,610,678` marked degree-`29` records.
The exact global remainder is therefore `42,581,669`, including every
five-or-more-black stratum.  Accordingly:

```text
four_black_status = PARTIAL
global_degree_29_status = OPEN
degree_29_witness_found = false
complete_four_black_no_go_claimed = false
complete_degree_29_no_go_claimed = false
```

## 5. Reproducibility

The exact producer and frozen archive are:

- [`certify_gram_belyi_m5_degree29_four_black_reduction.py`](../scripts/certify_gram_belyi_m5_degree29_four_black_reduction.py);
- [`gram_belyi_m5_degree29_four_black_reduction.json`](../data/enumerations/gram_belyi_m5_degree29_four_black_reduction.json).

One replay is casewise and bounded, launches no plane-tree enumeration, and
uses no timeout as mathematical evidence.  The next smallest residual target
is the `640` simple/simple `WW` master systems; the next largest uniform
target is the `1,920` `WW` systems with at least one critical endpoint.
