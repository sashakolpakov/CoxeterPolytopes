# Coupling the contact charge to the common height ratio

**Research proof candidate, pending independent review:** every finite-sided finite-volume real hyperbolic Coxeter polytope has dimension at most **166**. The new step is a local charge estimate of \(7/6\) when the common affine height ratio is at most \(21/20\). For larger ratios the frozen reviewed \(4/3\) estimate is already sufficient. No arithmetic hypothesis is used.

The [checker](../scripts/fresh_nikulin_coupled_charge_20261015.py) and [exact data](fresh_nikulin_coupled_charge_20261015.json) have passed generation and saved replay. This is producer verification, separate from independent proof review. All earlier [169 artifacts](fresh_nikulin_unit_charge_20261015_dimension169.md) and the [balanced unit-charge obstruction](fresh_nikulin_balanced_height_20261015.md) remain unchanged. There is no claim of publication priority, optimality, or formal proof-assistant verification.

## Local statement and common choices

Use the normalization and geometric hypotheses of the frozen [product-cusp proof](fresh_nikulin_narrow_20261015_products.md). Each affine factor has rank \(r\), Gram matrix \(G\) with root squares two, and positive kernel \(\delta\) normalized by \(\min\delta_k=1\). The prescribed omitted nodes are \(i,j\). The two nonnegative contact rows \(p,q\) are nonzero, and every nonzero entry is at least one. Each endpoint's total retained-support Schur cost is at most two. Set

\[
\sigma=\delta\cdot p,\quad \tau=\delta\cdot q,\quad
M=\max(\sigma/\tau,\tau/\sigma),\qquad
C=p_{S_i}^{T}G(S_i)^{-1}p_{S_i}
 +q_{S_j}^{T}G(S_j)^{-1}q_{S_j}.
\]

The normalized-row energy is \(E\ge0\), and

\[
W=\sigma\tau E/M=\min(\sigma,\tau)^2E.
\]

For the same preset omissions as in the [reviewed 4/3 lemma](fresh_nikulin_unit_charge_20261015_review.md), we prove

\[
\boxed{M\le21/20\quad\Longrightarrow\quad
 W+\tfrac76 C\ge T},
\qquad
T=\begin{cases}2r/9,&\text{ordinary factor},\\
2(r-1)/9,&\text{factor containing the forbidden narrow-facet node}.
\end{cases}
\tag{1}
\]

The choices are made before the terminating rows or their common ratio are known. In particular, ordinary E8 uses omissions `(1,6)`; the avoid-node E8 rules still use `(3,7)` or `(1,6)` as prescribed in the frozen proof. Both branches of the global argument use these identical choices.

This is compatible with the balanced counterexamples: at \(M=1\), the E8 common unit row has \(C=14/9\) and \(W=0\). A unit charge fails, while \((7/6)(14/9)=49/27>16/9\).

## The elementary height improvement

If the two rows have lower height bounds \(\sigma\ge h_1>0\), \(\tau\ge h_2>0\), then \(M\le B=21/20\) implies

\[
\min(\sigma,\tau)^2\ge
H(h_1,h_2):=
\max\!\left\{\min(h_1,h_2)^2,
\frac{\max(h_1,h_2)^2}{B^2}\right\}.
\tag{2}
\]

Indeed, the smaller height is at least the smaller individual lower bound, and is also at least the larger actual height divided by \(B\). No contact-label cutoff or approximation enters this estimate.

The case \(C\ge T/(7/6)\) proves (1) immediately. We therefore restrict to the strict small-cost range

\[
C<6T/7.
\tag{3}
\]

Every inverse retained Gram is entrywise nonnegative. A retained nonzero entry at node \(k\) costs at least its inverse diagonal \(c_{ik}\), and two entries in the same row also contribute their nonnegative cross term. These facts give the complete finite exceptional support classification below.

## Exceptional factors, including every forbidden node

The checker reconstructs the literal Grams, positive kernels, principal inverses and point resistances for A1, E6, E7, E8, F4 and G2. The six ordinary records and all 34 forbidden-node records exactly match the preset choices in the frozen proof.

For a single retained contact, orient the rows as

\[
p=a e_i+b e_l,\qquad q=c e_j,
\qquad b,c\ge1,\quad a=0\text{ or }a\ge1.
\]

Write \(R_{uv}\) for the energy of the normalized point-row difference \(e_u/\delta_u-e_v/\delta_v\). If \(a=0\), then \(E=R_{lj}\), while \(\sigma\ge\delta_l\), \(\tau\ge\delta_j\). Therefore

\[
W+\tfrac76C\ge
H(\delta_l,\delta_j)R_{lj}+\tfrac76c_{il}.
\tag{4}
\]

If \(a\ge1\), the heights are bounded below by \(\delta_i+\delta_l\) and \(\delta_j\). With

\[
u=\frac{\delta_l b}{\delta_i a+\delta_l b}\in[0,1],
\qquad
E(u)=(1-u)R_{ij}+uR_{lj}-u(1-u)R_{il},
\]

we obtain

\[
W+\tfrac76C\ge
H(\delta_i+\delta_l,\delta_j)\min_{0\le u\le1}E(u)
+\tfrac76c_{il}.
\tag{5}
\]

The energy is a convex quadratic. Its minimum is checked exactly at both endpoints and, if it lies in the interval, at its stationary point. All **54 single-contact occurrences** in the full 40-record list satisfy (4) and (5). The smallest recorded excess over the target is \(2363/1764>0\). Repeated forbidden-node choices are retained explicitly; the count is not a count of distinct diagrams.

The inverse-diagonal and cross-term census leaves exactly one two-contact pattern in (3): ordinary E7 with omissions `(3,6)`, one isolated A1 contact at node 4 in the first row and one at node 7 in the second. Its exact normalized energy is

\[
E(u,v)=\tfrac12(1+u^2+v^2)\ge\tfrac12.
\]

Both heights are at least one and the total cost is at least one, so

\[
W+\tfrac76C\ge\tfrac12+\tfrac76=\tfrac53>\tfrac{14}{9}.
\]

Every other two-contact support is excluded by its exact minimum cost, and three retained contacts exceed the threshold in every record. Thus this is complete support coverage with arbitrary real contact magnitudes satisfying the zero-or-at-least-one gap, not a scan of selected labels or unit rows. Zero retained contacts give the frozen point-row bound. The coincident-omission A1 avoid-node case has target zero and is immediate.

## Classical A: point surplus pays the excess linear loss

The following part is unconditional in \(M\). For affine A_r put \(m=r+1\), and choose opposite cycle nodes 0 and \(h=\lfloor m/2\rfloor\). Let

\[
R(k)=\frac{k(m-k)}m\quad(0\le k<m),\qquad
R_{ij}=R(h).
\]

For a retained node \(k\), its inverse diagonal is \(c_k=R(k)\), and its linear energy loss from the first omission is

\[
\ell_k=R_{ij}+R(k)-R(k-h),
\]

with arguments read modulo \(m\). Reflection of the cycle interchanges the omissions, so the same complete list covers the other row.

The frozen linear-loss identity, with nonnegative inverse entries and the contact gap, gives

\[
E\ge R_{ij}-\lambda C
\]

whenever \(\ell_k\le\lambda c_k\) at every possible retained contact. Since both heights are at least one, \(W\ge E\). In the range (3), only nodes with \(c_k<6T/7\) and \(c_k\le2\) can occur. Thus for \(\lambda\ge7/6\),

\[
W+\tfrac76 C\ge
R_{ij}-(\lambda-7/6)\frac{6T}{7}.
\tag{6}
\]

For each of the 18 cycle sizes \(m=2,\ldots,19\), the checker enumerates every node, takes the exact maximum of \(7/6\) and \(\ell_k/c_k\) over the allowed nodes, and verifies that (6) is at least \(T=2r/9\). This finite check covers all small ranks; it is not a cutoff on the theorem.

For larger ranks, the endpoint budget two implies that a retained contact is at cyclic distance at most two from its omission: for \(m\ge10\),

\[
\frac{3(m-3)}m>2,
\]

and the inverse diagonal increases up to the antipode. Write \(\varepsilon=m\bmod2\). At distance \(d\le2\),

\[
c=\frac{d(m-d)}m,\qquad
\ell\le\frac{(m+\varepsilon)d}m,
\qquad
6m\bigl(\tfrac76c-\ell\bigr)
\ge d(m-7d-6\varepsilon).
\]

The last expression is nonnegative for every even \(m\ge20\) and odd \(m\ge21\); all four distance/parity classes are checked symbolically. Finally the point surplus is positive:

\[
R_{ij}-2r/9=
\begin{cases}(m+8)/36,&m\text{ even},\\
(m-1)(m+9)/(36m),&m\text{ odd}.
\end{cases}
\]

This proves the ordinary estimate in every rank. Rotating the same omission pair to avoid a forbidden node preserves it and proves the weaker avoid-node target. A1 with coincident omissions is already handled separately.

## Classical B, C and D

The complete prefix-coordinate identities and their index-class coverage are inherited from the frozen [independent 4/3 review](fresh_nikulin_unit_charge_20261015_review.md). The new checker explicitly re-evaluates all **24** classes at threshold \(6T/7\).

For every non-distant contact the previously proved identity gives \(\ell\le c\), which is stronger than needed. The only other classes have costs \(r/4\), \(r/2\), \((r-1)/4\), or \((r-1)/2\), and remain strictly above the new small-cost threshold. For example,

\[
\frac r4-\frac{6}{7}\frac{2r}{9}=\frac{5r}{84}>0,
\]

and the modified B case is the same expression with \(r-1\). Each all-rank substitution is checked as an exact polynomial with nonnegative coefficients in independent nonnegative rank/index parameters. This is an inherited all-rank coordinate proof with new threshold comparisons, not an inference from a finite rank sample.

The point resistances \(r/4,r/2,(r-1)/4,(r-1)/2\) exceed their respective targets. Modified B2 is separate: its point energy is one, while a positive contact costs at least \(1/2\), whose \(7/6\) charge already exceeds the target \(2/9\). Consequently all B/C/D cases satisfy (1), without using the height-ratio condition.

## The global split at a common ratio

The frozen narrow-facet and actual-ridge arguments supply two terminating walls with a **single common ratio** \(M\) across all affine factors. The total endpoint costs satisfy \(\sum C_a\le4\). The sole avoid-node penalty remains \(2/9\), and \(N=\sum r_a=n-1\).

If \(M\le21/20\), summing (1) and applying the unchanged Lorentz identity gives

\[
-(x,y)\ge M\left(\frac N9-2\cdot\frac76-\frac19-2\right)
=M\frac{N-40}{9}.
\]

If \(M\ge21/20\), use the frozen 4/3 estimate with the same omissions:

\[
-(x,y)\ge M\frac{N-43}{9}.
\]

For \(N=166+t\), \(t\ge0\), the first lower bound is at least \(14+t/9\). The second is at least

\[
\frac{21}{20}\frac{123+t}{9}
=\frac{287}{20}+\frac{7t}{60}>14.
\]

Both contradict the strict narrow-facet bound \(-(x,y)<14\) for its norm-two framing normals. The ideal-vertex branch therefore excludes \(n\ge167\). The separately reviewed compact-facet branch gives \(n\le143\), so the combined candidate is \(n\le166\).

All geometric coverage conditions—finite-sided finite volume, cusp products, endpoint existence, common height ratio, actual framing ridges and the compact-facet alternative—are inherited from the reviewed 169 chain. The new proof is the smaller charge under a common-ratio condition and the two-branch optimization. Independent scrutiny of this new local step remains required.

```sh
python scripts/fresh_nikulin_coupled_charge_20261015.py
python scripts/fresh_nikulin_coupled_charge_20261015.py --verify
```
