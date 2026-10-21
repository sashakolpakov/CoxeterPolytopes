# Exact Wang correction crosswalk for the named project root lattices

**Status (2026-08-10): frozen exact positive-structure ledger.**  This adds
information that was absent from the earlier `51/18` yes/no ledger: the
correction lattice, a named Borcherds product or quasi-pullback, its weight,
the two divisor multiplicities when an `A1` component occurs, and twenty
same-form even-overlattice arrows between project rows.

The published inputs are Wang, Theorems 6.2 and 7.2 and Section 5, and
Gritsenko--Nikulin, Theorems 4.3, 6.1 and 6.2.  The project contribution is
the exact intersection with its named complete-root/divisor-lattice corpus.
The machine ledger contains 42 distinct named lattices: 26 corrected and 16
uncorrected.  Its rank-eleven--nineteen block is exhaustive and is exactly
`7+12=19`.

## 1. The missing exact high-rank rows

Write `S=U+L(-1)` and suppress minus signs in the table.  Except for the
explicit special rows below, the correction lives on

\[
                    M=U\oplus S=2U\oplus L(-1).
\]

`Borch(phi_0,L)` is Wang's pullback of the Jacobi input
`Theta_N/Delta`, equivalently the indicated quasi-pullback of Borcherds'
`Phi_12`.  Its principal multiplicity is `beta_0=1`.  If
`R(L)=mA1+R` and every component of `R` has Coxeter number `h`, then

\[
 \beta_1=2h-3,
 \qquad
 k=\left(12-\frac{\operatorname{rk}L}{2}-\frac{3m}{2}\right)h+12+3m.
\]

| `H` | hyperbolic lattice `S` | correction on `M` | weight; multiplicities |
|---:|---|---|---|
| 10 | `U+E8+A1` | `Borch(phi_0,E8+A1)` | `195`; `(beta_0,beta_1)=(1,57)` |
| 10 | `U+D8+A1` | `Borch(phi_0,D8+A1)` | `99`; `(1,25)` |
| 10 | `U+2D4+A1` | `Borch(phi_0,2D4+A1)` | `51`; `(1,9)` |
| 10 | `U+D4+5A1` | `Borch(phi_0,D4+5A1)` | `27`; `(1,9)` |
| 11 | `U+E8+A2` | **none** | Wang exception |
| 11 | `U+E8+2A1` | `Borch(phi_0,E8+2A1)` | `138`; `(1,57)` |
| 11 | `U+D8+2A1` | **none** | Wang exception |
| 11 | `U+2D4+2A1` | **none** | Wang exception |
| 12 | `U+E8+A3` | **none** | Wang exception |
| 12 | `U+E8+3A1` | **none** | Wang exception |
| 12 | `U+D8+3A1` | **none** | Wang exception |
| 13 | `U+E8+D4` | **none** | Wang exception |
| 13 | `U+D8+D4` | **none** | Wang exception |
| 13 | `U+E8+4A1` | **none** | Wang exception |
| 14 | `U+E8+D4+A1` | **none** | Wang exception |
| 15 | `U+E8+D6` | **none** | Wang exception |
| 16 | `U+E8+E7` | **none** | Wang exception |
| 17 | `U+2E8` | `Phi_{132,2E8}` | `132`; complete multiplicity-one divisor |
| 18 | `U+2E8+A1` | `Borch(phi_0,2E8+A1)` | `75`; `(1,57)` |

Thus the previous statement “seven corrected” is now an explicit list of
seven modular forms, not merely an existence flag.

## 2. Exact lower-rank project intersections

The frozen data also cross-links every complete-root row explicitly named in
the project's `H8` and `H9` ledgers and the four lower Bredon stock returns.
The ordinary corrected weights are:

```text
H8:
 U+E7          165       U+D6+A1       85
 U+D4+3A1       45       U+7A1         25
 U+A7           80       U+D7         114
 U+E6+A1        99

H9:
 U+E8          252       U+D8         124
 U+E7+A1       132       U+2D4         60
 U+D6+2A1       68       U+D4+4A1      36
 U+8A1          20

lower named stocks:
 U+D4       weights 72 and 40
 U(2)+D4   weights 40 and 24
 U+3A2     weight 39
```

The named exceptions in these rows are `U+D4+A2`, `U+D4+A3`,
`U+D5+A2`, and `U+E6+A2`.

There are two important non-split corrections which must not be silently
replaced by the ADE formula.

1. The project lattice `diag(2,-2^8)`, written
   `<2>+8A1(-1)`, is the sign reversal of
   `<-2>+8A1`.  Gritsenko--Nikulin's correction is on `U(2)+S` and is
   the additive lift

   \[
       \Delta_{5,D_7}=\operatorname{Lift}
       \left(\eta^3\prod_{i=1}^{7}\vartheta(z_i)\right),
   \]

   of weight `5`.

2. `U(2)+2D4` is exactly `U+N8`, where `N8` is the Nikulin
   index-two extension of `8A1`.  Hence
   `U+S ~= 2U+N8` and its correction is the complete-divisor
   quasi-pullback `Phi_{28,N8}` of weight `28`.

The dessin-ledger control now has its roles separated correctly:

```text
S = II_(1,9) = U+E8(-1)
actual Wang correction: U+S = 2U+E8(-1), Phi_{252,E8}, weight 252;

II_(1,9)+<2> ~= 2U+E7(-1), Phi_{165,E7}, weight 165
is instead the correction lattice for the H8 lattice U+E7(-1).
```

The two weight statements are both true, but they concern different
hyperbolic root systems.

Exact alternate project names are also frozen in the machine ledger:
`U+D10 ~= U+E8+2A1` (weight `138`),
`U+D6+3A1 ~= U+2D4+A1` (weight `51`),
`U+E7+2A1 ~= U+D8+A1` (weight `99`), and
`U+N8+A1 ~= U+D4+5A1` (weight `27`).

## 3. Same-form even-overlattice geometry

Every arrow below is an integral even inclusion of index two with the same
quadratic form; it is not a scaled-dual relation.

```text
H8:
 U+7A1 -> U+D4+3A1 -> U+D6+A1 -> U+E7
 U+A7 ---------------------------------> U+E7

H9:
                         U+D4+4A1
                        /          \
 U+8A1 ----------------              -> U+2D4 -> U+D8 -> U+E8
        \              /                         
         -> U+N8 = U(2)+2D4

 U+D6+2A1 -> U+E7+A1 -> U+E8

H10:
 U+D4+5A1 -> U+2D4+A1 -> U+D8+A1 -> U+E8+A1

H11:
 U+2D4+2A1 -> U+D8+2A1 -> U+E8+2A1

H12/H13:
 U+D8+3A1 -> U+E8+3A1
 U+E8+4A1 -> U+E8+D4 <- U+D8+D4
```

These arrows give two sharp cautions.

- Automorphic correction is **not downward-hereditary** under finite-index
  even sublattices: both lower nodes in the `H11` chain have no correction,
  while the maximal node has the weight-138 product.
- Arithmetic reflection maximality is **not** equivalent to automorphic
  correctability.  The companion exact computation proves that
  `U+E8+E7` is the maximal arithmetic reflection group in its rational
  class, while Wang proves that it has no automorphic correction.  Thus a
  proposed AG translation of maximality must retain information independent
  of the Borcherds divisor.

## 4. Scope boundary

Wang's theorem concerns the complete square-minus-two reflection group of
the lattice.  A finite-volume chamber cut by a chosen mixture of `(-2)`,
`(-4)`, divisibility-two, or scaled-dual mirrors is not thereby one of the
69 hyperbolic 2-reflective lattices.  In particular the project lattices
`U+E8(-2)+<-2>`, `U+E8(-2)+<-2>^2`, `U+A2(-1)^4`, and similar mixed-root
carriers are deliberately not assigned a Wang correction by this crosswalk.
The global complete-root theorem first tests whether their entire `(-2)`
arrangement has finite volume; only then is the `51/18` dichotomy applicable.

## 5. Reproducibility

- exact data:
  [`wang_project_complete_root_crosswalk.json`](../data/enumerations/wang_project_complete_root_crosswalk.json);
- lightweight certifier:
  [`certify_wang_project_complete_root_crosswalk.py`](../scripts/certify_wang_project_complete_root_crosswalk.py);
- focused tests:
  [`test_wang_project_complete_root_crosswalk.py`](../tests/test_wang_project_complete_root_crosswalk.py).

The certifier replays all Wang weight formulas, the `7/12` high-rank split,
the special non-split forms, and all twenty inclusion endpoints.  It performs
no root enumeration and stays far below the `0.5 GiB` cap.
