# Exact rank-three cusp directions with one omitted local facet

The remaining rank-three cusp calculation reduces to **91 exact local direction records**, including 44 represented over quartic real number fields. This is a complete local census for the stated case, not a global Coxeter-polytope construction or exclusion. The records count omission/target/direction combinations; no cross-field projective deduplication is claimed.

The literal source has three affine rank-three factors: cusp 4 with walls `[7,23,50,200]`, cusp 10 with `[7,63,119,197]`, and cusp 388 with `[32,39,162,181]`. A section that retains exactly three of a factor's four facets must have a bounded triangular Euclidean Coxeter link. Its irreducible type is affine A2, C2, or G2. There are respectively one, three, and six labeled angle patterns, so all four omissions in all three factors give 120 systems.

The [producer](../scripts/fresh_h20_rank3_omitted_sections_20261014.py) and [exact archive](fresh_h20_rank3_omitted_sections_20261014.json) enumerate these systems without an angle grid or floating-point predicates. There are 16 A2, 19 C2, and 56 G2 records. Their chosen real algebraic representations have degrees one, two, and four, with counts 33, 14, and 44. A representation's degree is not asserted to be the minimal projective field of its mirror.

## Conics give a complete direction classification

Let B be the positive-definite Gram matrix of the three retained old roots, and J its inverse. Choose planar outward root vectors u_i for the target Euclidean triangle. The implementation uses the following vectors, together with all distinct labeled angle patterns:

| Type | Planar vectors |
|---|---|
| A2 | `(2,0), (-1,sqrt(3)), (-1,-sqrt(3))` |
| C2 | `(1,0), (0,1), (-1,-1)` |
| G2 | `(1,0), (0,sqrt(3)), (-1,-sqrt(3))` |

The projected roots must be positive multiples of these rays, up to a common plane isometry. Set their relative lengths to `(1,x,y)`, with x,y strictly positive, and let X,Y be the two coordinate columns. Such columns represent an orthonormal plane in the old three-space, up to a common scale, exactly when

```
X^T J Y = 0,
X^T J X = Y^T J Y = q > 0.
```

After extracting the common square root from Y, the first two equations are rational conics in x,y. Conversely, those equations give a rank-two orthogonal projection, and

```
R = q B - X X^T - Y Y^T
```

is positive semidefinite of rank one. Choosing a nonzero column a of R recovers a normal with root pairings a. The script checks all rank-one minors, positive normal square, and the reconstruction identity exactly. Its kernel coefficients `k=J a` can be oriented strictly positive because the target triangle has a positive affine dependence.

For each system, the producer computes both rational resultants, removes their zero roots (which cannot have x,y positive), and adds their squarefree parts to the original conics. This preserves exactly their nonzero real solutions and removes multiplicities. A separating primitive element `t=x+c*y` then yields a rational univariate representation. The saved eliminants have degree at most four. Every real root of every irreducible factor receives a rational isolating interval; exact interval arithmetic decides the signs of x,y and all subsequent quantities. Substitution back into both conics is checked modulo the irreducible polynomial. Constant eliminants certify empty systems. No unresolved branch is silently discarded.

## The local parameter and the remaining work

For a record, let `e0=sum(k_i r_i)` and let `D=sum(m_i r_i)` be the source factor's positive affine relation. The normal is `e0+s D`; its square is independent of s. Normalize the horosphere by `D.z=-1`. On the retained inequalities, writing `w_i=-r_i.z>=0` gives `sum(k_i w_i)=-s`. The omitted inequality is redundant precisely when `sum(m_i w_i)<=1`, hence

```
-L <= s < 0,   L=min_i(k_i/m_i).
```

All polynomial normal coordinates, the normal square, positive kernel, null shift, L, and a minimizing wall are archived. The lower endpoint may create coincidences or additional incidence and requires separate treatment in a global construction.

The remaining task is to test the actual ridges of each full 210-inequality section over this interval. Passing the triangular cusp-link equations alone says nothing about those other ridges. The quartic directions in this census also show why the previous rational and quadratic direction banks could not claim completeness for this branch.

Both the primary saved-result replay and the [independent checker](fresh_h20_rank3_omitted_sections_20261014_review.md) pass. The latter imports no producer code: it derives all 120 original conics, directly saturates by `x*y` to remove zero coordinates, and checks the complete real solution sets and all 91 saved normal/parameter identities. It shares SymPy as its rational computer-algebra engine. The two methods agree on seventeen empty systems and every positive direction. Neither checks the global section ridges.

```sh
python3 scripts/fresh_h20_rank3_omitted_sections_20261014.py --verify
python3 scripts/fresh_h20_rank3_omitted_sections_20261014_review.py --verify
```
