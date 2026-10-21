# A subunit residual-cost envelope: dimension-164 candidate

**Research proof candidate, pending independent review:** the reviewed
[dimension-165 argument](fresh_nikulin_height_charge_20261015.md) sharpens to
dimension at most **164**. The new point is that two large E8 defects leave
less than one unit of cost for every other factor together. That residual
budget allows only one retained contact, whose total energy loss is at most
\(2/3\). This improves the aggregate loss from \(40/9\) to \(38/9\).
There is no new uniform unit-charge claim.

The [exact checker](../scripts/fresh_nikulin_defect_frontier_20261015.py) and
[certificate](fresh_nikulin_defect_frontier_20261015.json) bind the frozen 165
source and its [independent review](fresh_nikulin_height_charge_20261015_review.md).
They use the latter's standalone rational grounded-matrix helper read-only.
The new thresholds, complete low-cost singleton lists, classical tail
inequalities, and aggregate cases are recomputed. No earlier artifact is
modified. This package makes no publication-priority, current-best, or formal
verification claim.

## Unchanged local framework, with a larger height window

Use the same preset omissions, contact normalization, and geometric hypotheses
as the reviewed 165 proof. For each affine factor of rank r, put
\[
T=2r/9\quad\text{or}\quad2(r-1)/9
\]
in the unique factor avoiding the narrow-facet node. The energy is
\(W=\min(\sigma,\tau)^2E\), and C is the sum of the two retained-support
quadratic costs. The ratio
\(M=\max(\sigma/\tau,\tau/\sigma)\) is common across all factors.
Each endpoint's total cost is at most 2, hence \(\sum C\le4\).

The 165 exceptional support lists and small-cycle convex domains are
independent of the height-window endpoint. Re-evaluating all 71 singleton
energy bounds, the inexpensive two-row E7 case, and all 580 certified
small-cycle convex minima at
\[
\boxed{M\le B=21/20}
\]
still proves the same unit-charge classification. The A large-rank and B/C/D
proofs did not use height balance and remain unchanged. Thus every factor has
\(W+C\ge T\), except possibly the same 15 E7/E8 support occurrences:

| Type of actual failing factor | Target | Minimum cost | Defect upper bound |
|---|---:|---:|---:|
| Ordinary E8 | \(16/9\) | \(14/9\) | \(2/9\) |
| Ordinary E7 or avoid-node E8 | \(14/9\) | \(3/2\) | \(1/18\) |

An actual failing factor means \(T-W-C>0\). The table describes necessary
conditions, not an assertion that all listed support possibilities attain the
stated defects. At most two failing factors are possible because three would
consume at least \(9/2>4\) cost.

## The residual-product lemma

Under \(M\le21/20\), consider any collection of the prescribed affine
factors with total retained cost **strictly less than 1**. Then
\[
\boxed{\sum (T-W)\le2/3.} \tag{1}
\]
This includes ordinary and avoid-node factors of every rank.

First, for any positive definite retained Gram with root squares 2, its
inverse diagonal is
\[
J_{kk}=\frac1{2-v^TB^{-1}v}\ge\frac12.
\]
Its inverse is entrywise nonnegative. Every nonzero contact is at least 1, so
two retained coordinates, anywhere among the two rows of all the factors,
already cost at least 1. The strict residual budget therefore allows **at
most one retained coordinate globally**. A zero-cost factor has both rows
supported at their omitted nodes and satisfies the reviewed point-energy
bound \(W\ge T\). It remains to bound the loss of a single-contact factor.

### Exceptional factors

The checker reconstructs every retained inverse diagonal below 1 in all 40
ordinary/forbidden-node records. If \(T\le2/3\), inequality (1) for that
factor follows from \(W\ge0\). Otherwise the complete singleton
energy/height calculation gives the stronger \(W\ge T\) at \(B=21/20\).
Ordinary F4, whose target is \(8/9\), has no retained diagonal below 1.
There is no hidden multi-contact case because the total cost is strictly
below 1.

### All A cycles

For affine \(A_r\), set \(m=r+1\), choose omissions 0 and
\(h=\lfloor m/2\rfloor\), and write
\(R(d)=d(m-d)/m\). A retained node at distance d has inverse diagonal
\(R(d)\). Diagonal below 1 forces \(d=1\): for \(m\ge4\), every
\(2\le d\le m/2\) has
\(R(d)\ge2(m-2)/m\ge1\); m=2,3 already have only distance 1.

Write the single-contact row as \(p=a e_i+b e_k\), where k is a neighbor
of i, \(b\ge1\), and \(a=0\) or \(a\ge1\). The other row is supported
at its omitted node j.

If a=0, both heights are at least 1, and the worst neighbor is the one toward
the opposite omission. Hence
\[
W\ge R(h-1)=
\begin{cases}
m/4-1/m,&m\text{ even},\\
(m^2-9)/(4m),&m\text{ odd}.
\end{cases}
\]
Subtracting the ordinary target \(2(m-1)/9\) and adding \(4/9\) gives
\[
\frac{m^2+24m-36}{36m},\qquad
\frac{m^2+24m-81}{36m},
\]
respectively. These are nonnegative for all even m≥2 and odd m≥3, as verified
by their exact nonnegative-coefficient substitutions. Thus the loss is at
most \(4/9\), already smaller than the required \(2/3\).

If a≥1, one height is at least 2, so the smaller height is at least \(2/B\).
For m=2,3 the target is at most \(4/9\), so there is nothing further to prove.
For m≥4 the normalized-row energy is the squared distance from a point on
the segment between the two adjacent point potentials to the other omission's
potential. Its minimum, using the worse of the two adjacent directions, is
\[
E_{\min}=
\begin{cases}
\dfrac{m(m-2)}{4(m-1)},&m\text{ even},\\[3pt]
\dfrac{(m+1)(m-3)}{4(m-1)},&m\text{ odd}.
\end{cases}
\]
These follow by minimizing the exact convex quadratic
\(R(h)-\ell u+R(1)u^2\), with
\(\ell=(m+\varepsilon)/m\) and \(\varepsilon=m\bmod2\).
Both minima are at least \((m-3)/4\). Therefore
\[
W-T\ge\frac{m-3}{B^2}-\frac{2(m-1)}9
=\frac{302m-1102}{441}>0\qquad(m\ge4).
\]
Rotation proves the same ordinary estimate while avoiding a prescribed node
when m≥3. The coincident-omission avoid-A1 case has target zero. This covers
every A rank and every permitted omission case.

### All B/C/D ranks

The complete 24 prefix-coordinate classes are recomputed from the frozen
independent coordinate realization. Imposing inverse diagonal below 1 leaves
only the following possibilities:

| Family and position | Diagonal | Consequence |
|---|---:|---|
| Ordinary B, right omission, node 0 or 1 | \(r/4\) | r≤3, so T≤2/3 |
| Modified B, omission r−1, node 0 or 1 | \((r-1)/4\) | r≤4, so T≤2/3 |
| Modified B, omission r−1, node r | 1/2 | linear energy loss 0 |
| Separate modified B2 | at least 1/2 | T=2/9 |
| Modified C, omission 1, node 0 | 1/2 | linear energy loss 0 |
| Modified C, omission 1, node r | \((r-1)/2\) | r=2, so T=2/9 |

All other classes have diagonal at least 1. In each zero-linear-loss case the
nonnegative quadratic-remainder identity gives \(E\ge R_{ij}\), and both
heights are at least 1, so \(W\ge R_{ij}\ge T\). In the small-target cases,
\(W\ge0\) suffices. All modified choices are included via the previously
reviewed fork/path symmetries and the explicit modified B/C classes.

The constant \(2/3\) in (1) is sharp for this local domain. In ordinary
affine B3 with omissions (0,3), put both contact rows at node 0 with coefficient
1. Their normalized potentials agree, giving \(M=1\), \(E=W=0\), while
the sole retained cost is \(3/4<1\) and T=2/3. Its literal rationally scaled
Gram and inverse diagonal are checked exactly. This is a local witness, not
an assertion that an arbitrary product of extremal factors is realized by an
ambient polytope.

## Complete defect accounting

Let the total loss be \(L=\sum(T-W)\). If at most one factor fails unit
charge, then
\[
L\le\sum C+2/9\le38/9.
\]
If two fail, classify them by the number of ordinary E8 factors:

| Ordinary E8 among the two | Minimum combined cost | Remaining cost | Combined target | Total loss bound |
|---:|---:|---:|---:|---:|
| 0 | 3 | ≤1 | 28/9 | \(28/9+1=37/9\) |
| 1 | 55/18 | ≤17/18 | 30/9 | \(30/9+2/3=4\) |
| 2 | 28/9 | ≤8/9 | 32/9 | \(32/9+2/3=38/9\) |

The first row uses unit charge for all remaining factors, which do not fail
by definition. The other two use the strict-subunit residual-product lemma;
the individual failing factors themselves are bounded only by \(W\ge0\).
This proves the full envelope
\[
\boxed{L\le38/9\quad\text{whenever }M\le21/20.} \tag{2}
\]
It accounts for actual costs, not just a count of possible exceptional
support records.

## Dimension conclusion and scope

With \(N=n-1\), the sum of targets is \((2N-2)/9\), including the sole
avoid-node penalty. Equation (2) gives
\[
\sum W\ge\frac{2N-40}{9},\qquad
-(x,y)\ge M\frac{N-38}{9}\quad(M\le21/20).
\]
For \(M\ge21/20\), the unchanged four-thirds estimate with the same preset
omissions gives
\(-(x,y)\ge M(N-43)/9\). At \(N=164+t\), these are at least
\[
14+t/9,\qquad
14+7/60+7t/60,
\]
respectively. Both contradict the strict Nikulin narrow-facet pairing bound
\(-(x,y)<14\). Thus the ideal-vertex branch excludes n≥165; the inherited
compact-facet branch gives n≤143. Their combination is the candidate
\(\boxed{n\le164}\).

All geometric coverage—finite-sided finite volume, the narrow facet,
compact/ideal alternative, actual framing and terminating walls, common height
ratio, and endpoint Schur budgets—is inherited explicitly from the reviewed
165/169/166 chain. The new result is the widened local check and the complete
strict-subunit residual envelope. No arithmetic restriction or finite-label
cutoff is imposed.

For later work, a rank-residue argument alone does not rule out the relaxed
extremal combination: two E8 factors, one B3, one avoid-A1 and 24 zero-cost
ordinary E6 factors have total affine rank 164. Their individual endpoint cost
budgets may obstruct this combination, which is stronger information than
the summed cost used here. That is a separate future question, not part of
this proof or an achieved bound below 164.

```sh
python3 scripts/fresh_nikulin_defect_frontier_20261015.py
python3 scripts/fresh_nikulin_defect_frontier_20261015.py --verify
```

Generation and saved replay pass, as do source manifests, syntax, and local
links. Independent review of this new argument remains pending.
