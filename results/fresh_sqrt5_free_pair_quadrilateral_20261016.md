# Free pairing cannot close the chosen long edges into a quadrilateral

**All 210 pairs fail quadrilateral closure, even when the pairing of the two new walls is an arbitrary real number.** This concerns the specific uniform long-edge witnesses at the vertex D5 + 2E8. It does not exclude other edge parameters, longer faces, or compact dimension-21 polytopes.

The [exact checker](../scripts/fresh_sqrt5_free_pair_quadrilateral_20261016.py) and [certificate](fresh_sqrt5_free_pair_quadrilateral_20261016.json) retain the old 21-root Gram and both edges from the frozen [universal self-edge construction](fresh_sqrt5_universal_self_edges_20261016.md). Only the mutual product of the two duplicate roots is free. The result covers **144 pairs from different components** and **66 pairs within one component**.

## Continuous exclusion

Let `A` be the positive definite norm-two Gram of D5 + 2E8, put `C=A⁻¹`, and set

\[
T=2\phi^8,\qquad \phi^2=\phi+1.
\]

The new root at port `i` copies all old contacts except its contact with root `i`, which becomes `2−T`. For ports `i,j`, write `x` for the free new-to-new product. Eliminating **all 21 old roots** from the 23-root Gram gives

\[
\begin{pmatrix}p_i&x-b\\x-b&p_j\end{pmatrix},
\quad
p_i=2T-T^2C_{ii},\quad
p_j=2T-T^2C_{jj},\quad
b=A_{ij}+T^2C_{ij}.
\]

Both residual diagonals are strictly negative. There are 23 normals in an ambient vector space of dimension 22, so the residual must have rank at most one. A Lorentz realization also permits at most one negative direction.

Three intended spherical corners are copies of `A`. The fourth contains the two new norm-two roots. Even its two-root principal matrix can be positive definite only when `−2<x<2`; the acute Coxeter condition additionally requires `x≤0`. We relax these conditions to the entire closed real interval **[−2,0]**, without imposing any algebraicity or discrete angle condition.

The residual determinant is

\[
D(x)=p_ip_j-(x-b)^2.
\]

For every pair, the checker proves exactly that `D(−2)>0` and `D(0)>0`. These **420 endpoint checks** certify the whole interval, because

\[
D(-2+2t)=(1-t)D(-2)+tD(0)+4t(1-t)>0
\quad(0\le t\le1).
\]

Thus the residual is negative definite throughout the finite-corner interval. The full matrix has **rank 23 and signature (21,2)**. Every possible fourth spherical corner is excluded before conjugate positive semidefiniteness, integrality or Coxeter-label restrictions could help. The fixed individual edges retain their positive conjugate residual diagonals; no omitted conjugate case can repair this physical signature failure.

The smallest of the 420 exact endpoint determinants is

\[
20169647+32635176\phi>0,
\]

at ports 2 and 3 with `x=−2`. The proof is a continuous interval certificate, not a sample of candidate labels.

## What the rank equation still permits

The two real solutions of the rank condition are

\[
x_\pm=b\pm\sqrt{p_ip_j}.
\]

The same exact inequalities show `x₋<−2` and `x₊>0`. The positive solution violates the acute wall convention. The negative solution makes the new walls **ultraparallel** rather than incident at a fourth vertex.

The negative branch does retain the original vertex and its two finite neighbors as a real partial configuration. This is checked separately, so failure of the quadrilateral is not attributed to a missing endpoint. In the orthogonal decomposition with old roots spanning the positive space, let `t₀²=−1`, let `omega_i` be the old Gram-dual vector, and put `a_i=sqrt(−p_i)>0`. Use

\[
\text{new}_i=e_i-T\omega_i+a_it_0,\qquad
v_i=t_0-\frac{a_i}{TC_{ii}-1}\omega_i.
\]

Then the original vertex is `t₀`, the neighbor `v_i` has square `−1/(TCii−1)²`, and its active diagram is again D5 + 2E8. For both directions of each pair, the checker verifies

\[
(-p_j)(TC_{ii}-1)^2
  -(-p_i)T^2C_{ij}^2>0.
\]

These **420 strict sector checks** prove that `v_i` lies strictly inside the other new wall. All three vertices lie in the same time cone and both specified finite edges remain present.

This last construction uses declared real square roots. It does not certify arithmeticity, a common integral lattice, a full edge star, or finite volume. A two-face completion preserving these two edges would need additional side walls and hence at least five sides. The calculation does **not** rule out such a completion.

Generation and saved replay passed, including **54 recursive manifests and 526 source references**. The source, certificate and report are frozen.

```sh
python3 scripts/fresh_sqrt5_free_pair_quadrilateral_20261016.py
python3 scripts/fresh_sqrt5_free_pair_quadrilateral_20261016.py --verify
```
