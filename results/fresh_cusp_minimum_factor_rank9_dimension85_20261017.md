# Dimension at most85 when every cusp factor has rank at least9

**Theorem.** Every finite-sided, finite-volume real hyperbolic Coxeter
polytope whose irreducible affine cusp factors all have rank at least9
has dimension at most85. There is no arithmetic hypothesis or limit on
the number of factors at a cusp. Compact polytopes satisfy the cusp
condition vacuously. A polytope with a smaller cusp factor is outside
this theorem; the unrestricted ceiling remains163.

The proof combines the [distance-seven corner payments](fresh_cusp_factor_payment_framework_20261017.md),
the [exact cusp-link curvature inequality](fresh_cusp_link_shape_curvature_20261017.md),
and the [existing complete internally checked163 proof](fresh_nikulin_two_budget_frontier_20261015_dimension163.md).
The finite arithmetic step is a complete integer partition enumeration,
not a finite-rank extrapolation.

## Exact local comparison

At a cusp with factor ranks r_i, put d=sum r_i=n−1,
H(t)=product_i(1+t+...+t^r_i)=sum h_j t^j, V=H(1), and
G=sum h_j min(j,d−j). The curvature is

` kappa(L)=(5/2)G−dV/2. `

The minimum-rank9 hypothesis excludes affine subdiagrams of order at most9.
The inherited small-core proof therefore pays every triangle at least2
and every quadrilateral at least1 using graph distance at most7.
At a finite vertex the cost is at most7n.

In an affine factor of order r+1, the number of paid pairs is at most
7(r+1), and trivially at most r(r+1)/2. Each pair occurs in exactly
V/(r+1) link-edge fibers. Thus the complete cusp cost is at most

` V * sum_i min(7,r_i/2). `

The [exact enumerator](../scripts/fresh_cusp_minimum_factor_rank9_dimension85_20261017.cpp)
checks every nondecreasing integer partition of each d from85 through162
with parts at least9. In all50,330,195 such partitions it proves

` kappa(L)/V − sum_i min(7,r_i/2) > 1/5. `

All arithmetic is integral. If c=sum_i min(14,r_i), the condition is
exactly

` 5*(5G−(d+c)V)−2V > 0. `

The complete enumeration includes50,378,399 nonempty partitions with
total rank at most162; smaller totals are traversal prefixes and are
not asserted to satisfy the inequality. In particular the six-factor
profile(14,14,14,14,14,14) in ambient dimension85 has a negative margin
for this coarse payment comparison. It is only a local model, not an
asserted hyperbolic chamber.

## Exhaustive dimension argument

Suppose n>=86. The inherited unrestricted theorem first gives n<=163.
For every integer n in this finite remaining range, beta_n>7n, where
beta_n=n(n−6)/8 for even n and (n−1)(n−5)/8 for odd n.
The curvature inequality gives

` D5 >= beta_n v_f + sum_v kappa(L_v) + 5(n−1). `

The corner payments give

` D5 <= 7n v_f + sum_v V_v sum_i min(7,r_vi/2). `

The exact local comparison makes these inequalities contradictory.
There is at least one vertex; the positive additive constant supplies
additional slack. This excludes every n>=86 and proves the theorem.
Use of the existing163 proof is explicit: no assertion about an
unexamined rank tail is made by the finite enumeration.

## Arithmetic completeness and verification

The recursive generator appends a rank r at least as large as the last
rank and visits every prefix. Therefore it generates each unordered
rank partition exactly once. Appending r multiplies H by1+t+...+t^r;
a sliding-window integer convolution computes every new coefficient.
At each partition, symmetry and sum h_j=product(r_i+1) are checked.

For every r>=9 in the tested range, (r+1)^2<=2^r. Consequently V<=2^81,
every h_j<=V, and G<=81V. All intermediate magnitudes in the enumerator
are below2^94, safely within the signed128 integer type. No floating
point calculation participates in the certificate. The compiler command
keeps assertions enabled.

The [Python wrapper](../scripts/fresh_cusp_minimum_factor_rank9_dimension85_20261017.py)
compiles and runs the exact enumerator and independently counts every
rank partition by the unbounded-coin recurrence. It also checks all
finite-vertex margins and the exact dimension85 local obstruction.
The [saved certificate](fresh_cusp_minimum_factor_rank9_dimension85_20261017.json)
binds the sources and inherited proof packages. A saved replay recompiles
and reruns the enumeration; it does not merely reload the saved counts.

The independent source/composition review is recorded in the
[review note](fresh_cusp_minimum_factor_rank9_dimension85_20261017_review.md).
Generation and a fresh saved replay both terminate with PASS, rerunning all
50,378,399 exact convolutions. The independent review and its separate
Ferrers-conjugation partition counter also pass, matching every rank total;
it does not claim a duplicate convolution implementation. Source syntax,
direct bindings and local links pass. The source, data and report are final.
External mathematician refereeing, proof-assistant formalization and
publication priority have not been established. No new polytope is
constructed by this theorem.
