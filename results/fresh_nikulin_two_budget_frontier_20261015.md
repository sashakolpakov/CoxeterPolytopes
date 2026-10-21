# Ordered endpoint budgets: dimension-163 candidate

**Internal research proof candidate, pending independent review:** the frozen,
independently reviewed [164 argument](fresh_nikulin_defect_frontier_20261015.md)
sharpens to dimension at most **163**. The two endpoint Schur budgets are each
at most two. Keeping their order, rather than replacing them by their sum,
gives the complete near-balanced energy-loss envelope
\[
\boxed{\sum_a(T_a-W_a)\le4\quad\text{if }M\le21/20.}\tag{1}
\]
The proof does not assert a uniform unit-charge estimate. Such estimates have
known balanced counterexamples. Nor does it extrapolate from a finite rank or
finite dihedral-label search.

The [checker](../scripts/fresh_nikulin_two_budget_frontier_20261015.py) and
[certificate](fresh_nikulin_two_budget_frontier_20261015.json) preserve all
previous files. They use the frozen standalone Fraction grounded-matrix
[review helper](../scripts/fresh_nikulin_height_charge_20261015_review.py)
read-only. New support lists, convex minima, radical bounds, tail inequalities,
and ordered packing decisions are recomputed. Published geometric inputs and
the compact-facet branch remain the explicitly inherited reviewed chain.
No publication-priority, current-best, or formal-verification claim is made.

## Framework and the residual lemma

Use exactly the preset omissions of the reviewed 164/165/169 packages. For a
rank-r affine factor, let
\[
T=2r/9,\qquad T=2(r-1)/9\text{ in the sole avoid-node factor}.
\]
Its endpoint contact rows p,q have nonzero entries at least one. The positive
affine kernel \(\delta\) is normalized by \(\min\delta_i=1\). Put
\[
\sigma=\delta\cdot p,\quad \tau=\delta\cdot q,\quad
M=\max(\sigma/\tau,\tau/\sigma),\quad
W=\min(\sigma,\tau)^2E.
\]
Here E is the squared distance of the normalized Euclidean point-potential
means. The same M applies to every factor. Write \(C_0,C_1\) for its two
retained-support quadratic costs. The actual endpoint geometry gives
\[
\sum_a C_{a,0}\le2,\qquad\sum_a C_{a,1}\le2.\tag{2}
\]

The new residual-product lemma is:

| Residual endpoint budgets | Maximum total loss \(\sum(T-W)\) |
|---|---:|
| \((10/9,4/3)\) | \(20/9\) |
| \((1,3/2)\) | \(22/9\) |
| \((3/2,1)\) | \(22/9\) |

It covers arbitrary numbers and ranks of prescribed factors, with at most one
avoid-node factor, for every \(M\le21/20\). All three maxima in the relaxed
local domain are attained by the exact table records below. The proof first
establishes a complete continuous envelope for every factor with
\(C_0,C_1\le3/2\), then performs finite exact ordered-budget packing.

## Complete continuous domains at finite ranks

For an omitted node i and a retained node k, the inverse diagonal is
\(\delta_k^2R(i,k)\), where R is the squared distance of the normalized
point potentials. The retained inverse Gram is entrywise nonnegative.
Therefore every actual support lies among all subsets of nodes with diagonal
at most \(3/2\), and its cost is bounded below by evaluating every nonzero
coordinate at one, including all same-row cross terms. No support-size cutoff
is used.

For every pair of retained supports, the checker includes both possibilities
for each omitted entry: zero, or at least one. Empty full rows are discarded
because their positive affine height must be nonzero. If both retained
supports are empty, the previously reviewed point-energy inequality gives
\(W\ge T\), so they contribute no positive loss.

For full supports U,V, the normalized means lie in the two finite convex
hulls. The checker constructs an exact rational closest-pair certificate.
For their signed Euclidean Gram K, barycentric vector x and the two simplex
constraints A, it verifies
\[
x\ge0,\quad Ax=(1,1),\quad
Kx+A^t\lambda\ge0,\quad
x_i(Kx+A^t\lambda)_i=0.
\]
The explicit Euclidean models make K positive semidefinite. These identities
prove the global minimum \(E\ge e=x^tKx=-\lambda_0-\lambda_1\), including
singular minimizer sets. There is no floating-point optimization in the
certificate.

The full supports also give height lower bounds
\(s=\sum_{k\in U}\delta_k\), \(t=\sum_{k\in V}\delta_k\). Thus
\[
W\ge H_B(s,t)e,\qquad
H_B(s,t)=\max\{\min(s,t)^2,\max(s,t)^2/B^2\},\quad B=21/20.
\tag{3}
\]
Square roots in s,t are either exact rational numbers or replaced by certified
64-bit rational lower bounds. H is increasing in each argument, so this is a
valid downward energy bound. Surviving support-cost floors are checked to be
exactly rational; radical lower bounds used while excluding other supports
cannot discard a feasible support.

The finite coverage is A cycles of sizes 2 through 17, B/C ranks 2 through 7,
D ranks 4 through 17, all ordinary/avoid E6/E7/E8/F4/G2 records, and the
separate modified B2. This gives **91 models, 1,116 retained-support pairs and
3,630 continuous full-support domains**. Exactly 1,395 computed loss upper
bounds are positive. Each is dominated, with the same modified/ordinary flag,
by one of the following 28 rational records: the dominating record has no
larger cost floor in either endpoint and no smaller loss upper bound.

| Kind | Cost floors \((C_0,C_1)\) | Loss upper bound |
|---|---|---:|
| ordinary A1 | \((1/2,0)\) | \(2/9\) |
| ordinary A2 | \((2/3,0)\) | \(4/9\) |
| ordinary A4 | \((4/5,4/5)\) | \(8/9\) |
| ordinary B2, B3, B4, B5 | \((0,r/4)\), r=2,3,4,5 | \(2r/9\) |
| ordinary D4, D5 | \((r/4,0)\), r=4,5 | \(2r/9\) |
| ordinary E6 | \((0,4/3)\), \((4/3,0)\) | \(4/3\) |
| ordinary E7 | \((1/2,1)\), \((1,1/2)\) | \(14/9\) |
| ordinary E8 | \((8/9,2/3)\) | \(16/9\) |
| modified B3 through B7 | \(((r-1)/4,0)\), r=3,…,7 | \(2(r-1)/9\) |
| modified E6 | \((0,4/3)\) | \(10/9\) |
| modified E7 | \((0,3/2)\), \((2/3,5/6)\), \((5/6,2/3)\) | \(4/3\) |
| modified E8 | \((1,1/2)\), \((8/9,2/3)\) | \(14/9\) |
| modified F4 | \((0,1)\), \((1/2,3/4)\) | \(2/3\) |
| modified G2 | \((0,1/2)\) | \(2/9\) |

The surviving table happens to consist of identical point rows with W=0.
That observation is **not** the completeness premise: every continuous
support domain is separately bounded by (3) and its saved convex certificate.
The JSON records the dominating table index for every positive-loss domain.

## Classical tails, with no rank cutoff

The following models are the same verified prefix-coordinate realizations as
in the frozen unit-charge proofs. Let \(s_k=e_1+\cdots+e_k\).
For B, use \(v_0=0,v_1=e_1,v_k=s_k/2\) for k≥2; the kernel squares are
\((1,1,4,\ldots,4,2)\). For C, use \(v_k=s_k\), metric one half the
usual Euclidean metric, and kernel squares \((1,2,\ldots,2,1)\).
For D use the B prefix vectors through r−2 and
\(v_{r-1}=(s_{r-1}-e_r)/2,v_r=s_r/2\), with kernel squares
\((1,1,4,\ldots,4,1,1)\).

Inserting these formulas into \(\delta_k^2R(i,k)\le3/2\) gives exactly
the retained-node lists below. Nodes not listed have larger diagonal, so
positive cross terms cannot restore them. The omission itself is optionally
present in each full support.

| Type and ranks | Omissions | Possible retained nodes, in order | Fixed-coordinate bound for E | Height-square bound |
|---|---|---|---:|---:|
| B, r≥8 | (0,r) | {1}; {r−1} | \((r-2)/4\) | \(2/B^2\) |
| modified B, r≥8 | (r−1,0) | {r−2,r}; {1} | \((r-3)/4\) | \(2/B^2\) |
| C, r≥8 | (0,r) | {1}; {r−1} | \((r-2)/2\) | 1 |
| modified C, r≥8 | (1,r) | {0,2}; {r−1} | \((r-3)/2\) | 1 |
| D, r≥18 | (0,r) | {1}; {r−1} | \((r-2)/4\) | 1 |

For example, every convex mean in the left B support has zero coordinates
2 through r−1, while every right mean has those coordinates equal to one
half. This gives \((r-2)/4\). Its right full support has minimum kernel
entry \(\sqrt2\), so the common ratio gives the stated height bound.
The other rows follow from the same fixed coordinates. Exact affine
polynomials in r−8 or r−18 show that W−T is nonnegative in every row; only
the D18 starting value is zero. These are all-rank identities, not tests at a
large numerical cutoff. Literal boundary cases are also checked in the JSON.

For the remaining A cycles set m=r+1, \(h=\lfloor m/2\rfloor\),
\(\varepsilon=m\bmod2\). For m≥18 even or m≥19 odd, a diagonal
\(d(m-d)/m\le3/2\) forces d=1. Two retained coordinates in one row
already cost more than \(3/2\). Thus each full support consists of its
omission and at most one adjacent node.

If both full supports are singletons, their cyclic distance is at least h−2,
and \(W\ge R(h-2)\). Its excess over the ordinary target has numerator
\[
m^2+8m-9(4+\varepsilon)^2
\]
over the positive denominator 36m. If either full support has two nodes, one
height is at least two. The frozen exact linear-loss identity gives
\(E\ge R(h)-2(m+\varepsilon)/m\), hence
\[
W-T\ge
\frac{(9-2B^2)m^2+(2B^2-72)m-81\varepsilon}{9B^2m}.
\]
Substituting m=18+2t or 19+2t gives positive constant terms and nonnegative
remaining coefficients in both numerators. This proves **W≥T**, a stronger
statement than the old tail's W+C≥T; the two are not conflated.

Rotating opposite A omissions avoids any prescribed node when m≥3. Fork
swaps do the same in D. Ordinary B/C omissions already avoid all other nodes;
their exceptional choices are the modified rows above. In these cases using
the stronger ordinary target and treating the record as ordinary only
enlarges the residual envelope. Modified A1 has target zero and needs no
positive-loss record. All small boundary cases are in the finite census.

## Exact packing and the actual failing-factor split

The table has minimum positive total cost at least one half. An exact finite
recursion takes nondecreasing multisets of its records, respecting both
ordered budgets and allowing at most one modified record. It recomputes the
three residual maxima in the first table in 110, 151 and 163 memoized states.
The witnesses are respectively an ordinary E8 plus B2, an ordinary E7 plus
two B2, and an ordinary E7 plus ordinary B2 and modified B3. This is a
complete finite packing calculation, not a linear relaxation or sampled set.

Now call a factor actually failing if \(T-W-C_0-C_1>0\). The frozen 164
proof, already valid at B=21/20, restricts such factors to its 15 listed
E7/E8 support occurrences. Two are E7 same-row supports at nodes4 and7.
The other endpoint row must then be supported on its omission. Exact convex
energy and height bounds, with and without the same-row omitted entry, give
W+3/2>14/9 in all four cases. These two occurrences cannot actually fail.

The checker reconstructs all remaining 13 rows and their ordered floors:

| A selected actually failing factor | Ordered cost floors | Target | Residual budgets | Residual loss |
|---|---|---:|---|---:|
| ordinary E8 | \((8/9,2/3)\) | \(16/9\) | \((10/9,4/3)\) | \(20/9\) |
| ordinary E7 or modified E8 | \((1/2,1)\) or \((1,1/2)\) | \(14/9\) | \((3/2,1)\) or \((1,3/2)\) | \(22/9\) |

If there is no actually failing factor, summing unit charge gives total loss
at most four. Otherwise select **any one** such factor. Bound its own loss by
its target, using W≥0, and apply the complete residual lemma to every other
factor. Each row gives total loss at most four. This argument also includes
additional actually failing factors in the residual product. If the selected
factor was modified, still allowing one modified residual record is a harmless
enlargement. Equation (1) follows.

## Dimension consequence and sharp local accounting

The targets sum to \((2N-2)/9\), where N=n−1 and the subtraction is the sole
avoid-node penalty. Therefore (1) gives
\[
\sum W\ge(2N-38)/9,\qquad -(x,y)\ge M(N-37)/9
\quad(M\le21/20).
\]
For M≥21/20, the unchanged four-thirds estimate with the **same preset
omissions** gives \(-(x,y)\ge M(N-43)/9\). At N=163+t the two lower bounds
are respectively at least
\[
14+t/9,\qquad14+7t/60.
\]
They contradict the strict narrow-facet inequality \(-(x,y)<14\) for every
t≥0. The ideal-vertex branch therefore gives n≤163; the inherited compact
facet branch gives n≤143. Combined, the research candidate is
\(\boxed{n\le163}\), without arithmetic hypotheses.

The local loss constant four is sharp. Take two ordinary E7 factors with both
rows at node7, and two ordinary B2 factors with both rows at node0. Their
common M is one, W=0, the endpoint cost totals are exactly (2,2), and the
total target is four. Adding a zero-target modified A1 and 24 zero-cost E6
factors with W=T gives total rank163. This is a checked configuration of the
local charge domain, **not** an ambient polytope. Further improvement may
need rank-sensitive accounting or additional geometric restrictions; no
claim of an optimal dimension bound follows from this local witness.

Generation and saved replay both **PASS**; syntax, source manifests and local
links also pass. The primary source, certificate and this report are frozen.
The
published narrowness, actual framing/terminating-wall existence, common
height ratio, endpoint budgets and compact-facet coverage remain inherited
from the frozen reviewed chain; this certificate addresses the new residual
envelope and global accounting.
