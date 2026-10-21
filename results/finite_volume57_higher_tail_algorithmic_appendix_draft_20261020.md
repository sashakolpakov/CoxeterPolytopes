# Higher-dimensional exclusions: payment families and exact algorithms

This is an insertion draft for the existing finite-volume paper, not another
paper or a new computation. It replaces a chronological list of 27 reports
by the common mathematical argument, its bounded domains, and the separate
all-rank argument. All arithmetic quoted below is from existing saved runs;
no producer, checker, import body, or build was executed in preparing it.

The scalar definitions and rules from dimension66 onward are written out,
with the finite preferred-image tables declared as exact mathematical
supplement data in Section8. Dimensions59–65 are being recovered separately
from this ambient66 base. Aggregate costs never define a marked rule.
The remaining geometric proof dependencies are identified explicitly.

## 1. One comparison proves all bounded intervals

Let the ambient dimension be n and put d=n−1. Write D_5 for the sum of
5 minus the number of sides over all two-faces. A simultaneous corner
payment rule w assigns nonnegative weights to actual corners, pays at
least2 on a triangle and at least1 on a quadrilateral, and therefore gives

    D_5 ≤ U_n v_f + Σ_(ideal v) C_w(v).                         (1)

Here U_n is an upper bound for the cost of every rank-n spherical forest.
The geometric curvature inequality in the body of the paper gives

    D_5 ≥ β_n v_f + Σ_(ideal v) κ(v) + E_n,  E_n>0,             (2)

where

    β_n = n(n−6)/8                         if n is even,
          (n−1)(n−5)/8                    if n is odd.         (3)

The weaker positive correction E_n≥5(n−1) suffices in this appendix; the
stronger two-extrema correction used in the low-dimensional certificates
is compatible with it. Thus the following is the only comparison needed.

**Tail comparison.** If one and the same rule has the stated local payments,
U_n≤β_n, and κ(v)−C_w(v)>0 at every ideal vertex, dimension n is impossible.
Indeed subtraction gives

    0 ≥ (β_n−U_n)v_f + Σ_(ideal v)(κ(v)−C_w(v)) + E_n >0.

This also covers a polytope with no ideal vertices. A saved finite maximum
or cusp table alone does not establish the local-payment premise of(1).

For an ideal vertex let r_1,…,r_s be its affine factor ranks, so Σr_i=d.
Set

    H(t)=∏_i(1+t+⋯+t^r_i)=Σ_(j=0)^d h_j t^j,
    V=∏_i(r_i+1),   G=Σ_j h_j min(j,d−j),
    κ=(5G−dV)/2.                                             (4)

If a connected affine factor of rank r has at most (r+1)c(r) paid marked
pairs, the actual factor-fiber incidence theorem gives

    C_w(v) ≤ V Σ_i c(r_i).                                   (5)

The normalized surplus to be checked is consequently

    M_c(r_1,…,r_s)=(5G−dV)/(2V)−Σ_i c(r_i).                  (6)

Formula(5), including the multiplicity V/(r+1), comes from the geometric
factor-fiber theorem. It is not inferred just from the graph pair count.

## 2. Exact terminal enumeration and why it is complete

The following describes a finite algorithm, including its mathematical
domain. It applies to all the bounded intervals in Section4.

For each fixed n and envelope c, a terminal profile is a nondecreasing
list of positive ranks summing to d and satisfying one of the terminal
conditions below. Enumerate it recursively: choose the next rank at least
the previous one and at most the remaining rank; recurse with the
remaining rank decreased by that choice. Reject an extension only when
its multiplicity violates an explicitly stated terminal condition. Accept
exactly when the remaining rank is zero. This enumerates each terminal
profile once, and halts because every choice decreases a nonnegative
integer by at least one.

For an accepted profile compute H by integer convolution, starting at
H=1 and convolving successively with the all-one list of length r+1.
Compute V independently as the product of these lengths. Check symmetry,
degree d, and Σh_j=V. Compute G both from its defining weighted sum and,
independently, from cumulative probabilities of the integer sum. Evaluate
(6) with exact fractions. A certificate contains every terminal profile
once and the exact values of H,V,G,Σc,M_c. A verifier reconstructs the
domain, checks set equality and uniqueness, and then checks these
identities. An upper-bound-only certificate need not prove attainment;
any claimed minimum must additionally carry an attaining profile.

The domain reduction uses the concentration move

    (a,b) → (a−1,b+1),   1≤a≤b,                              (7)

discarding a zero part. The centered sum of uniforms on these two rank
intervals increases in convex order. Equivalently its finitely supported
stop-loss sums increase; the exact two-factor kernel proof is given in
the preceding cusp-algorithm appendix. Since

    κ/V=(3d−5 E|2X−d|)/4,

the normalized curvature does not increase. Apply(7) only where

    c(a−1)+c(b+1) ≥ c(a)+c(b).                               (8)

Then M_c does not increase either. Also
(a+1)(b+1)−a(b+2)=b−a+1>0, so the link vertex count decreases.
Every move increases Σr_i² by 2(b−a+1); as Σr_i²≤d², no sequence is
infinite. This proves termination of the reduction independently of the
enumeration. Positive terminal normalized surplus implies positive
surplus for every original cusp. If a terminal absolute minimum is used,
positivity is checked first, so multiplication by the larger original V
has the correct direction. Negative terminal absolute margins cannot be
extended in this way.

Three useful terminal forms are:

* **Low-L only:** at most one part in1,…,L−1, all other parts at leastL.
  No concentration of high parts is allowed. This is essential for the
  varying high-rank envelopes used at n59–62.
* **Linear/cap at L:** arbitrary copies of L, at most one part belowL,
  and at most one part aboveL. Concentrate pairs both belowL or both
  aboveL. Here c(r)=min(r/2,L/2).
* **Piecewise:** list a set R of freely repeated ranks and a set O of
  ranks allowed at most once, together with at most one rank below the
  first repeated rank and one above the last. The permitted moves are
  (7) in the linear low segment and constant high segment, and the
  separately verified duplicate-rank moves displayed in Section4.

For the finite-vertex bound U_n, use the spherical-type alphabet and the
ordinary/contextual forest-cost dynamic program of the preceding
algorithm appendix. Every admitted type has positive integer rank. Each
state records remaining rank and all multiplicities or guard bits that
affect a future contextual charge. Its Bellman certificate checks every
last-item upper inequality. An equality predecessor gives attainment when
needed. Enumerating all fixed large-component multiplicities before this
recurrence covers interactions among them and the residual types. Neither
one representative forest nor an unproved compression of the state
constitutes this finite-vertex bound.

## 3. Two intrinsically defined high-tail rules

A feature is a node of degree at least3 or an endpoint of a nonsimple
bond, including the parallel bond. Distances are unweighted graph
distances. Conditions on a whole connected component refer to the whole
actual endpoint diagram, not a selected subdiagram.

### 3.1 The feature-eight rule, n98–108

At a finite vertex pay a unit to pairs of distance at most9 and to every
pair in a whole connected component of order at most17. The finite
forest bound is U_n=9n. At an ideal vertex pay a unit for distance at
most8, for whole component order at most17, or for distance9 with at
least one mark within14 steps of a feature. These conditions are united,
not added. In an affine factor of order m≥18, the exact graph totals are

|Family|Total paid pairs|
|---|---:|
|A cycle|8m|
|C path|min(9m−45,8m−4)|
|B, one fork and a doubled far end|min(9m−44,8m−2)|
|D, two forks|min(9m−43,8m)|

For example on a C spine there are 8m−36 pairs at distance at most8;
the extra distance9 starting positions occupy the first16 and last16
places, giving min(m−9,32). The B and D spines give respectively first15
and last16, or first15 and last15, positions; their short leaves add2
or4 pairs. This proves the formulas at every rank. At order at most17,
even all pairs cost at most8m. Hence c(r)=min(8,r/2).

The local replacement theorem must be stated separately: each changed
actual corner either has the required retained short/feature path, or
has a preferred whole component of order at most17. This is the
distance-eight overlap/routing theorem, not a consequence of the graph
count just proved. Its complete core and arbitrary-white proof is an
internal geometric input to this family.

### 3.2 Distance k, n109–163

For k=9,…,13, pay a unit at both finite and ideal vertices for distance
at mostk or whole component order at most17. The ideal count is at most
min(km,binom(m,2)) for orderm; thus c(r)=min(k,r/2). The finite bound is
U_n=kn. The families are used on

    k=9:109–116;  k=10:117–126; k=11:127–136;
    k=12:137–146; k=13:147–156.

For n157–163 use k=14 and the distance-only rule; U_n=14n and
c(r)=min(14,r/2). The small-component clause is unnecessary there.
The connected finite/affine classification gives the ideal pair bound
by counting distances on a spine and its at most four short leaves.
The actual-face payment is supplied by the complete distance-k local
replacement theorem. It includes connected cores, overlapping cores,
and arbitrary supporting whites; it is not a statement that every
two-face is visible in one fixed small diagram.

In these six cases use the linear/cap terminal domain with L=2k.
Existing exact calculations give the following positive minima over all
terminal profiles and all ambient dimensions in the indicated interval:

|n|Envelope cap|Minimum normalized surplus|
|---|---:|---:|
|98–108|8|112101/1419857|
|109–116|9|45873/2476099|
|117–126|10|145234/4084101|
|127–136|11|40603/559682|
|137–146|12|172096/1953125|
|147–156|13|830407/9565938|
|157–163|14|1476616/20511149|

For every row β_n−U_n>0. These scalar results, together with the local
theorems just identified, exclude all n98–163 by Section1.

## 4. Lower tail: envelope families, not a chronology of files

This section supplies the scalar part of n59–97. Section8 supplies the
marked rules for n66–97; the ambient59–65 layer ledger is a separate
pending insertion.
All values below are normalized per affine-factor order as in(5).

The basic ideal-five rule pays distance at most5 and distances6–9 with
one mark within5 steps of a feature. Above order17 its all-rank totals are

    p_A(m)=5m,
    p_C(m)=5m−15+Σ_(k=6)^9 min(m−k,14),
    p_B(m)=5m−6 +Σ_(k=8)^11 min(m−k,13),
    p_D(m)=5m+3 +Σ_(k=10)^13 min(m−k,12).                    (9)

Define g(r)=max_X p_X(r+1)/(r+1). The spine-and-leaf proof counts
starting positions in the union of initial and final intervals; it is
valid for every order, including overlapping intervals. The slopes of
the tree numerators are successively9,8,7,6,5. Their endpoints show

    g(r)≤163/23,  with equality at r=22,
    g(r)=(5r+56)/(r+1) for r≥24.                             (10)

The low-only profiles at n59–62 use g(r) itself above rank16. They do
not replace it by163/23. The latter replacement would invalidate the
saved positivity check in this range.

For n63–80 the ideal envelopes are piecewise-linear at low ranks and
constant at163/23 from rank17 onward. The exceptional values are:

|n|Low linear segment|Values after that segment, through rank16|
|---|---|---|
|63|r≤9|54/11,21/4,147/26,6,94/15,107/16,118/17|
|64|r≤9|109/22,16/3,74/13,6,94/15,107/16,118/17|
|65|r≤10|16/3,74/13,6,94/15,107/16,7|
|66|r≤10|16/3,75/13,6,32/5,107/16,120/17|
|67–69|r≤10|65/12,75/13,43/7,32/5,107/16,120/17|
|70|r≤11|75/13,43/7,32/5,107/16,120/17|
|71|r≤11|75/13,87/14,97/15,107/16,120/17|
|72–74|r≤11|75/13,87/14,97/15,111/16,120/17|
|75|r≤12|87/14,97/15,111/16,126/17|
|76–77|r≤12|87/14,20/3,111/16,126/17|
|78–80|r≤12|87/14,20/3,111/16,127/17|

The entries in the last column are in consecutive ranks, beginning
immediately after the stated linear segment, where c(r)=r/2.

The terminal domains are specified without any search cutoff:

|n|Freely repeated ranks R|Additional exceptional ranks O|
|---|---|---|
|59–62|all ranks≥9|at most one rank1,…,8; no high concentration|
|63–64|9,…,17|at most one rank1,…,8 and one rank≥18|
|65|10,…,17|at most one rank1,…,9 and one rank≥18|
|66|10,11,12,14,15,16,17|13, plus one low<10 and one high>17|
|67–69|10,11,13,14,15,16,17|12, plus one low<10 and one high>17|
|70–74|11,13,14,15,16,17|12, plus one low<11 and one high>17|
|75–80|12,14,16,17|13,15, plus one low<12 and one high>17|

Each rank in O is allowed at most once. The special duplicate moves
2j→(j−1,j+1) have normalized-cost gains

    n66:j13,11/65; n67–69:j12,23/1092;
    n70:j12,19/182; n71–74:j12,16/91;
    n75:j13,4/105 and j15,7/2040;
    n76–77:j13,5/21 and j15,83/408;
    n78–80:j13,5/21 and j15,107/408.

For example the rank12 duplicate move is not valid at n66: its cost
gain is negative. This is why rank12 remains freely repeated in that
row. The stated domains avoid this invalid reduction.

Three final envelope shapes suffice for n81–97:

|n|Envelope c(r)|Terminal repeated/exceptional ranks|
|---|---|---|
|81–86|r/2 to13; c14=20/3,c15=111/16,c16=127/17; c≥17=95/13|R={13,16,17}; O={14,15}, one low<13, one high>17|
|87–91|r/2 to14; c15=7,c16=127/17; c≥17=217/29|R={14,16,17}; O={15}, one low<14, one high>17|
|92–97|min(15/2,r/2)|linear/cap at15|

The extra duplicate gains for the first row are5/48 at rank14 and
107/408 at rank15; for the second,8/17 at rank15. Thus all are covered
by Section2's terminating reduction.

The finite upper bounds and saved strictly positive scalar outputs group
as follows. A minimum displayed over an interval is the minimum over
all dimensions and all terminal profiles in that interval.

|n|Finite upper U_n|Minimum normalized ideal surplus|
|---|---|---:|
|59|783/2|1/340|
|60|404|31/6000|
|61|420|3/250|
|62|434|2741/371712|
|63|449|3505/201344|
|64|464|2385/146432|
|65|479|61/1936|
|66|495|13367/495924|
|67–68|504,512 respectively|851/165888|
|69|544|3641/13608|
|70|8n|1819/51840|
|71|8n|23641/338688|
|72–74|8n|6875/49152|
|75|8n|63/74698|
|76–77|8n|11615/448188|
|78–80|9n|409611/2171546|
|81–86|9n|940675/11294304|
|87–91|9n|146417/911250|
|92–97|9n|10121/4194304|

Every finite bound in this table is at mostβ_n, including the equality
cases. The first seven rows additionally have absolute ideal surplus
at least 91/2,54,64,73,139/46,240/23,445/23, respectively. Those stronger
numbers are not needed when using(1)–(2).

## 5. The all-rank end: geometric reduction to two ordered budgets

The remaining n≥164 case does not follow from a finite list of cusp
partitions. It uses the published narrow-facet theorem: a finite-volume
Coxeter polytope has a facet whose actual framing walls, with norm-two
normals, satisfy the strict bound

    −(x,y)<14.                                               (11)

The exact framing/endpoint argument is part of the geometric spine:
terminating walls can be chosen at finite or ideal endpoints of actual
outgoing edges. At a cusp of this facet write N=n−1 for the sum of
affine ranks. Each factor has two omissions chosen before the contact
rows. In the sole factor containing the narrow-facet node the omissions
avoid it, and the target is reduced by2/9. No other factor is modified.

In a factor of rankr let δ be its positive affine kernel, normalized
by minδ=1. For the two actual endpoint rows p,q write positive null
heights σ=δ·p and τ=δ·q. Their ratio q_0=σ/τ is the same in every
factor. Put M=max(q_0,q_0^−1). The two Schur costs obey independent
budgets

    Σ C_0≤2,  Σ C_1≤2.                                     (12)

These are ordered budgets; replacing them by their sum loses the sharp
residual estimate. For normalized Euclidean point potentials define
E to be the squared distance of the two normalized row potentials and
W=min(σ,τ)^2 E. The Lorentz identity is

    −(x,y)= (M/2)ΣW −q_0−q_0^−1,
    Σ targets T = (2N−2)/9,
    T=2r/9 ordinary, T=2(r−1)/9 in the modified factor.      (13)

The geometric proof supplies common height and actual framing walls;
neither is asserted just from the local matrix computations below.

### 5.1 The far-height range

For the prescribed omissions, at every affine rank,

    W+(4/3)(C_0+C_1) ≥ T.                                  (14)

Here is a reproducible proof specification for this internal local lemma.
Let R_ij=||v_i−v_j||². The point-potential expansion gives

 E=R_ij−Σw_k(R_ij+R_ik−R_kj)
         −Σz_k(R_ij+R_jk−R_ki)+||a−b||².                  (15)

Weights w,z are convex, each at most1, and a nonzero contact magnitude
is at least1. The retained inverse has nonnegative entries, so its
cost dominates all diagonal unit-contact costs. For B,C,D use the
coordinate potentials in Section6; split a row into its endpoint,
prefix, and fork classes. In each class either

    R_ij+R_ik−R_kj ≤ (4/3)δ_k²R_ik,

or δ_k²R_ik>3T/4, which by itself pays the target. The finite class
split has24 rows, but its coverage is the endpoint/prefix/fork partition,
not an assertion about24 sampled ranks. In each prefix row substitute
k=k_0+u,r=r_0+u+v with u,v≥0; the cleared difference has nonnegative
coefficients. For A, R(d)=d(m−d)/m. Below the trivial cost threshold,
small cycles have at most one adjacent contact; for m≥10, the endpoint
budget≤2 excludes distances≥3. The remaining cleared difference is

    d(m−4d−3ε)/(3m),  ε=m mod2,

for the admissible d=1,2. It is nonnegative. The finitely many A1 and
exceptional-type omission cases use all subsets of retained nodes below
the trivial cost threshold; convex quadratic minimization covers their
arbitrary magnitudes, including an omitted-node contact either0 or≥1.
The only cheap two-contact support has energy(1+u²+v²)/2. This gives
(14) with no finite-angle cutoff.

For precision, the classical24-class partition in this check is the
following. A contact range is a separate symbolic class; isolated listed
nodes are separate classes. Empty ranges contribute nothing.

|Type and omission row(i,j)|Retained contact classes|
|---|---|
|B ordinary(0,r)|1; 2≤k≤r−1; r|
|B ordinary(r,0)|0; 1; 2≤k≤r−1|
|B modified(r−1,0)|0; 1; r; 2≤k≤r−2|
|B modified(0,r−1)|1; r; 2≤k≤r−1|
|C ordinary(0,r)|1≤k≤r−1; r|
|C modified(1,r)|0; 2≤k≤r−1; r|
|C modified(r,1)|0; 1≤k≤r−1|
|D ordinary(0,r)|1; 2≤k≤r−2; r−1; r|

The other C ordinary row follows by path reversal; the other D row and
avoid-node choices follow by fork automorphisms. In each interval write
k=k_min+u and r=k+v+ε with the fixed endpoint gapε, u,v≥0.
Compute R directly from(19), cost=δ_k²R(i,k), and
loss=R(i,j)+R(i,k)−R(k,j). The exact certificate checks either the
nonnegative coefficients of(4/3)cost−loss after clearing its positive
denominator, or those of cost−3T/4 with a positive constant term.
This gives an explicit finite symbolic algorithm and its full index
domain. Modified B2, where prefix endpoints coincide, uses the separate
three-node matrix in Section6. For the exceptional cases enumerate
all cheap supports in the grounded matrices of Section6 with the
prescribed ordinary/avoid-node omissions; the low-target modified A1
case has T=0 and is immediate. No rank is omitted by these reductions.

Consequently for M≥21/20, (12)–(14), together with
q_0+q_0^−1≤2M, give

    −(x,y) ≥ M(N−43)/9.                                    (16)

The exact finite24-row identities and the40 exceptional omission rows
are internal certificate inputs; their independent coordinate derivation
is in the four-thirds review identified in Section9. Omitting these
identities from a final paper would leave a local lemma to be proved,
not a published theorem.

### 5.2 The near-height range

For 1≤M≤B=21/20 the sharper ordered-budget assertion is

    Σ(T−W)≤4.                                              (17)

Sections6–7 define its entire finite and infinite domain and the exact
verification algorithm. With(13), it implies

    −(x,y)≥M(N−37)/9.                                      (18)

If N=163+t,t≥0, (18) is at least14+t/9. In the far range(16) is at
least14+7t/60. Both contradict(11). Thus the cusp branch excludes
every n≥164.

If the chosen narrow facet has no cusp, it is compact. The separate
compact-facet payment proof gives ambient n≤143. One coarse all-rank
check is to compare the finite curvature in dimensiond=n−1 with
17(d+1). For d=143+2t their difference is t²/2+36t+3/2; for
d=144+2t it is t²/2+73t/2+19. Both are positive. The compact-facet
payment theorem, not the cusp payment rule localized without proof,
supplies this upper bound. These exhaustive branches finish the
unbounded end.

## 6. Exact continuous-support domain in the near range

This is the finite-model algorithm behind(17). It has no contact-label
cutoff and includes all arbitrary magnitudes at least1.

Its endpoint cap3/2 applies only to the residual product after selecting
a failing factor in Section7. The original product still has budgets2;
the selected factor is classified separately without this residual cap.

For classical rankr set s_k=e_1+⋯+e_k. The potentials and kernel squares
are

    B: v0=0,v1=e1,vk=s_k/2(k≥2),
       δ²=(1,1,4,…,4,2);
    C: vk=s_k, metric one half the ordinary metric,
       δ²=(1,2,…,2,1);
    D: v0=0,v1=e1,vk=s_k/2(2≤k≤r−2),
       v_(r−1)=(s_(r−1)−e_r)/2, v_r=s_r/2,
       δ²=(1,1,4,…,4,1,1).                                 (19)

Write R(i,j)=||v_i−v_j||². For a cycle of orderm,
R(i,j)=d(m−d)/m, d=(i−j) modm, and δ=1. Ordinary omissions are
(0,floor(m/2)) for A and(0,r) for B,C,D. Modified B uses(r−1,0),
modified C uses(1,r). A rotation or D fork automorphism avoids a
prescribed node without changing the ordinary geometric estimate;
weakening its target by2/9 is safe.

The finite list consists of A ranks1–16; B ranks2–7, with modified
B at3–7; C ranks2–7, both ordinary and modified; D ranks4–17;
E6,E7,E8,F4,G2 ordinary and each avoid-node choice; and one modified
B2 case. Exceptional node labels are fixed as follows. A star has
central node0 and consecutive labels along its arms, in the listed
order. Use arm lengths(2,2,2),(1,3,3),(1,2,5) for affine E6,E7,E8.
Their integral symmetric matrices have diagonal2 and edge−1. For F4,
use diagonal(2,2,2,4,4) and consecutive edges(−1,−1,−2,−2), with
scales(1,1,1,2,2); for G2 use

    [[2,−1,0],[−1,2,−3],[0,−3,6]], scales(1,1,3).

A matrix A with these scales has physical Gram
A_ij/sqrt(s_i s_j). If k>0 spans kerA, normalize
δ_i²=s_i k_i²/min_j(s_j k_j²). Ground node0 and solve A u=f for
f=e_i/k_i−e_j/k_j, u_0=0; then R(i,j)=min(s_l k_l²) f·u.
These equations uniquely determine the rational resistance values.

Ordinary exceptional omissions are E6:(2,4), E7:(3,6), E8:(1,6),
F4:(0,4), G2:(0,2). For forbidden nodef choose:

* E6: the first two of(2,4,6) different fromf;
* E7:(2,5) if f∈{4,7}, otherwise(4,7);
* E8:(1,6) if f∈{3,7}, otherwise(3,7);
* F4:(0,4) if f∈{1,3}, otherwise(1,3);
* G2: the two nodes other thanf.

Modified B2 uses matrix
[[2,0,−2],[0,2,−2],[−2,−2,4]], scales(1,1,2), and omissions(0,1).
Its physical forbidden node is2. A saved producer metadata field equal
to0 in this one model is a placeholder; it must not be interpreted as
the forbidden node when stating the avoid-node theorem.

For each omissiono, eligible retained nodes satisfy δ_k²R(o,k)≤3/2.
Enumerate every subsetS of those nodes and compute its minimum cost

    c_o(S)=Σ_(k∈S) δ_k²R(o,k)
       +Σ_(k<l∈S) δ_kδ_l[R(o,k)+R(o,l)−R(k,l)].             (20)

Discard a subset only if this lower bound exceeds3/2. All cross terms
are nonnegative. The actual retained inverse cost for any magnitudes
at least1 is therefore at least(20). For every pair of surviving
subsets, independently include or omit the omitted-node contact on
each side; included means magnitude at least1. Require each full row
to be nonempty. The two-empty-retained case has the direct zero-cost
bound W≥T and need not enter the residual list.

For full supportsU,V form their two convex hulls of potentials. Their
squared distance has a positive semidefinite signed Gram matrixK and
two constraints Ax=(1,1),x≥0. The algorithm enumerates active subsets
and solves their rational KKT equations. A returned witness is checked
by

    Ax=(1,1), x≥0, Kx+A^tλ≥0,
    x_i(Kx+A^tλ)_i=0, e=x^tKx=−λ_1−λ_2.                   (21)

These identities prove that e is the global minimum over the entire
continuous domain: expanding at x gives a nonnegative quadratic and
nonnegative linear slack. Singular active systems may be skipped;
a finite model is accepted only once some witness satisfying(21) is
provided. Exhausting the finite active-set list always halts; failure
to find a witness is a failed certificate, never silent pruning.

Let s=Σ_(U)δ_k and t=Σ_(V)δ_k. Since the height ratio is at mostB,

    W≥max(min(s,t)²,max(s,t)²/B²) e.                         (22)

Actual unbounded magnitudes are covered by normalized convex weights,
and their positive lower bounds give s,t. Thus(21)–(22) certify a loss
upper bound T−W on a full support domain. The computations use rational
outward lower bounds for square roots in(20),(22). Whenever a support
survives, its cost is also checked exactly; a radical lower bound is
not silently substituted for an exact surviving cost.

The complete saved calculation has91 finite models,1116 retained
support pairs,3630 nonempty continuous domains, and1395 positive-loss
domains. These counts are checksums. Completeness follows from the
model list, classification, every-subset construction, both omitted
contact states, and the all-rank tails below.

For the remaining ranks the diagonal criterion reduces the supports
directly. The following are simultaneous lower bounds for every point
of both full convex hulls:

|Type|Ranks|Omissions|Retained supports on the two sides|Energy lower bound|Height-square lower bound|
|---|---|---|---|---:|---:|
|B ordinary|r≥8|(0,r)|{1};{r−1}|(r−2)/4|2/B²|
|B modified|r≥8|(r−1,0)|{r−2,r};{1}|(r−3)/4|2/B²|
|C ordinary|r≥8|(0,r)|{1};{r−1}|(r−2)/2|1|
|C modified|r≥8|(1,r)|{0,2};{r−1}|(r−3)/2|1|
|D ordinary|r≥18|(0,r)|{1};{r−1}|(r−2)/4|1|

Sets in this table are allowed retained nodes, not a requirement that
all be contacted. A row can additionally contact its omitted node.
Every difference of full-support points has a fixed middle block of
coordinates of magnitude1/2 in B,D or1 in C. Counting those coordinates
proves the energy bound. Its product with the height bound minusT
is affine in the excess rank, nonnegative at the initial rank and
strictly increasing. This is an all-rank proof, not three boundary
tests extrapolated to larger ranks.

For A, m≥18 even or m≥19 odd, cost≤3/2 permits only one adjacent
retained node per side. If both full supports are singletons, their
cyclic separation is at least floor(m/2)−2, and W−T has lower numerator

    m²+8m−9(4+ε)² over36m.

If either full row has two nodes its height is at least2/B; the lower
numerator is

    (9−2B²)m²+(2B²−72)m−81ε over9B²m.

Substitute m=18+2t or19+2t, t≥0. Both polynomials have nonnegative
coefficients and positive constant term. Therefore all these tails have
W≥T and contribute no residual loss.

## 7. Ordered packing and the loss-four theorem

For every positive-loss finite domain replace its two cost floors and
loss upper bound by a dominating record: smaller or equal cost floors,
larger or equal loss, and the same modified flag. The exact Pareto list
has the following28 records; a range denotes one record for each rank.

|Flag/type|(C_0,C_1) floors|Loss upper bound|
|---|---|---:|
|ordinary A1|(1/2,0)|2/9|
|ordinary A2|(2/3,0)|4/9|
|ordinary A4|(4/5,4/5)|8/9|
|ordinary B_r,r2–5|(0,r/4)|2r/9|
|ordinary D_r,r4–5|(r/4,0)|2r/9|
|ordinary E6|(0,4/3),(4/3,0)|4/3 each|
|ordinary E7|(1/2,1),(1,1/2)|14/9 each|
|ordinary E8|(8/9,2/3)|16/9|
|modified B_r,r3–7|((r−1)/4,0)|2(r−1)/9|
|modified E6|(0,4/3)|10/9|
|modified E7|(0,3/2),(2/3,5/6),(5/6,2/3)|4/3 each|
|modified E8|(1,1/2),(8/9,2/3)|14/9 each|
|modified F4|(0,1),(1/2,3/4)|2/3 each|
|modified G2|(0,1/2)|2/9|

The dominating map is checked for every one of the1395 positive-loss
domains; a favorable list of records without this map is insufficient.
The table is not assumed realizable by one polytope.

Define F(i,p,q,a) as the maximum total loss of a nondecreasing multiset
of records with indices at leasti, ordered cost budgets(p,q), and at
mosta modified records. Its recurrence is the maximum of0 and

    loss_j+F(j,p−C_0j,q−C_1j,a−flag_j)                       (23)

over all admissible j≥i. Each record has C_0+C_1≥1/2. Thus a recursion
has length at most floor(2(p+q)); all data are rational with a common
denominator, so the state set is finite. Induction on this remaining
length proves that(23) covers every multiset and that the checked
upper/equality-predecessor certificate gives its exact maximum.
The ordered results are

    F(0,10/9,4/3,1)=20/9,
    F(0,1,3/2,1)=F(0,3/2,1,1)=22/9.                        (24)

The choice of a failing factor uses a second complete local result:
if T−W>C_0+C_1, then after two explicitly excluded same-row E7
supports, it has one of the following cost floors and targets:

    ordinary E8: (8/9,2/3), T=16/9;
    ordinary E7 or modified E8: (1/2,1) or(1,1/2), T=14/9.

Its finite proof enumerates every cheap retained support of the same
grounded exceptional matrices. The two E7 supports at nodes4,7 of
one row have direct convex certificates, including both omitted-node
states, proving W+3/2>T; they are excluded before this conclusion.
The remaining13 cases have precisely the displayed floors. This
failure classification is another required certificate, distinct from
the28-record residual domination.

If no factor fails, total loss≤Σ(C_0+C_1)≤4. Otherwise choose one
failing factor. In the E8 case, the remaining budgets are at most
(10/9,4/3), so its own loss at most16/9 plus(24) is at most4.
In the E7-type case the remaining budgets are(3/2,1) or(1,3/2),
so14/9+22/9=4. All additional failures are already included in that
residual maximization. Allowing one modified record in the residual
even when the chosen factor is modified only weakens the bound.
This proves(17).

## 8. Complete marked-rule ledger for n66–97

This ledger is a mathematical definition of the rules. A finite supplement
table used below is a list of unordered node pairs, not an executable
module or an aggregate total. Its schema and reconstruction are specified
in8.6. All node labels in this section start at0.

### 8.1 Labelled targets and elementary rules

For a classical affine diagram of orderm, use these labelled edges:

* A: the m-cycle, all labels3;
* C: path0,…,m−1, with label4 on its first and last bonds;
* B: spine0,…,m−3 with its last bond labelled4, and simple leaves
  m−2,m−1 attached to0;
* D: spine0,…,m−5 with simple leaves m−4,m−3 attached to0 and
  m−2,m−1 attached to m−5.

Every unlisted off-diagonal entry is zero. The finite A_r path is
0,…,r−1; finite B_r has label4 only on its last bond; finite D_r
has spine0,…,r−3 and leavesr−2,r−1 attached tor−3.
The affine A1 parallel pair and the exceptional affine E/F/G diagrams
of order at most9 are wholly paid in these families.

Let q_b be the ideal basic rule with unit weight for distance≤b, or
distance b+1,…,9 with a mark withinb steps of a feature. At small
orders add the precise preferred-image tables below, rather than
paying all pairs automatically. The union operation assigns weight1,
not one unit per witness. Define the D cross-fork set at orderm by

    X_m={unordered pairs at distance≥m−4}
          ∪{the pair of the two degree-three nodes}.         (25)

For the relevant m≥14 this has nine elements. The modified small
tables already remove exactly these nine pairs when that removal is
specified below.

Define finite rules R9 and R8:

    R9: unit on same-component distance≤9 or whole rank≤17;
    R8: unit on same-component distance≤8, or distance9 with
        a whole-component leaf mark, or whole rank≤17.        (26)

All other same-component weights are zero until an override is stated.
Across different components the base weight is zero. All the role
conditions below refer to whole finite components.

### 8.2 Two finite-nine rules cover n78–97

Use R9 and add a half-unit between any mark of a component of rank≤4
and a mark of either D12,…,D16 at its fork center or a short leaf,
or A13/A15 at its unique middle. This is the complete broad rule;
count a marked pair once. Its total is at most9n for n≤108.

Use this rule for n81–97. It slightly enlarges the earlier n87–97 rules, which used
only selected A1/A2 cross roles. Enlargement preserves all local
payments, and the independently proved broad budget remains9n.
Thus the same finite rule can be used throughout n81–97, without
transferring an unproved earlier equality or changing any ideal cost.

For n78–80 use instead the wider finite-eight rule F8 of8.3, including
its D11 and A17 coupons, and enlarge its same-component part from R8
to R9. Its cost is at most9n: F8 costs at most8n for n≤80, and the
only added pairs in a rankr≥18 component are nonleaf distance9 pairs,
of which there are r−11 in A/B and r−12 in D, each at mostr. Ranks≤17
already pay all pairs. Summing this increment over all components
costs at mostn. The wide coupons cannot be dropped in this interval:
the ideal-five cycle replacements can require the A17 and D11 roles.

The corresponding ideal tables are:

|n|Whole small orders|Modified small orders|All larger orders|
|---|---|---|---|
|92–97|m≤16: all pairs|m17: P7(17,X)|q7|
|87–91|m≤15: all pairs|m16: P7(16,X); m17:P7(17,X)|q7|
|81–86|m≤14: all pairs|m15,16,17:P6(m,X)|q6|
|78–80|m≤13: all pairs|m14,15,16,17:P5(m,X)|q5|

Here X is A,B,C,D and the finite pair sets P_b are the exact preferred
unions defined in8.6, including their validated cycle and fork removals.
Their totals are respectively

    P7(17):119,111,108,127;
    P7(16):112,107,103,111;
    P6(15):90,100,95,96;
    P6(16):96,107,103,111;
    P6(17):102,111,108,127;
    P5(14):70,87,84,82; P5(15):75,100,95,96;
    P5(16):80,107,103,111; P5(17):85,111,108,127.             (27)

The totals corroborate the pair tables; they do not define them.

### 8.3 Base finite rule at n76–77 and its enlargements

Use R8. Add a half-unit between any mark in a component of rank≤5 and
a mark of D11,…,D16 at the fork center or short leaf, of A13/A15 at
the middle, or of A17 at path position5 or11. Also add a half-unit
between an A4 endpoint and the unique long-arm neighbor of the D13
fork center. These role sets are disjoint where their types coincide.
Call this rule F8. Its complete finite budget is8n for n≤80.

At n70–75 use F8 with two further coupons: upgrade the half-unit for
A1/D12 with a D12 short-leaf mark to1, and add a half-unit for either
A2 mark paired with the short terminal node of B10 incident to the
doubled bond. Its budget is8n for n≤75. Choosing this enlarged rule
also at n75 is safe, although that particular historical proof did
not need both additions.

At n69 make three additional changes: delete the finite A17 endpoint
pair(0,16); upgrade A1/D11 short-leaf pairs from half to1; and give
weights1 and1/2, respectively, to A1/A11 and A2/A11 with the A11
mark its middle node5. Its exact finite maximum at n69 is544.

At n67–68 replace the same-component rank17 tables of A,B,D by the
finite sets T17(X) defined in8.6; retain all the preceding coupons.
Their costs are114,105,127. A17 remains coupon-bearing; B17/D17
are not identified with it. The finite maxima are504 and512.

At n66, retain this rule, change the short-leaf distance9 pairs in
whole D_r,r≥20, from1 to1/2, and add a half-unit between any mark
of a whole component of rank≤5 and a short leaf of whole D10.
Finally upgrade the existing A2/D11 short-leaf half-unit to1.
The exact finite maximum is495. These are actual marked rules;
the representative finite-cost algorithm may allow B17 the larger
D17 cost127 while its actual cost is105.

The last rule has the following convenient doubled forest upper-cost
formula, suitable as a starting state for the n59–65 layer extraction.
Let b(t) be the doubled base and j(t,u) the doubled component-pair cost.
For ranks≤16, b=r(r−1); at17 use A228 and B/D254 as upper values;
for r≥18 use16r−68 for A/B,16r−64 for D, and16r−66 for D≥20.
The actual B17 base is210. For each orientation(small,large), add to j:

    2 rank(small), if rank(small)≤5 and large=D10;
    3 rank(small), if rank(small)≤5 and large=D11,…,D16;
      rank(small), if rank(small)≤5 and large=A13 orA15;
    2 rank(small), if rank(small)≤5 and large=A17;
    2, for (A1,D11),(A1,D12),(A2,B10),(A4,D13),
           (A1,A11),(A2,A11);
    4, for (A2,D11).                                        (28)

The extra2 or4 encode upgrades or distinct roles, not repeated
charging of the same pair. The literal role definitions above determine
the actual weights. In the counting proof every other type is replaced
by a same-rank representative only after checking base and every
cross-cost domination. In particular, order17 has two distinct states.

### 8.4 Exact ideal-zero and half-weight patches below n78

Start from P5 at orders14–17, all pairs at orders≤13, and q5 at
orders≥18. Apply all patches whose row is at least the current
ambient dimension. A patch sets a weight; it does not subtract its
value from an already zero entry.

|Used for n|Patch|
|---|---|
|66–77|Set affine D,order17,pair(1,11) to0.|
|66–75|Set affine B,order15,pairs(0,11),(11,13),(11,14) to0.|
|66–74|At order13 use the four explicit tables in8.5; at D,order17,set the set H17 below to1/2.|
|66–71|At D,order16,set H16 below to1/2.|
|66–70|Set B,order14,pair(0,10) to0 and B,order15,pair(1,11) to0.|
|66–69|At order12 use the four explicit tables in8.5.|
|66 only|At B,order12,set(9,10),(9,11) to0; at C,order12,set(0,10),(1,11) to0; at B,order14,set(10,12),(10,13) to0.|

The fixed half-weight sets are

    H17={(0,10),(2,12),(0,11),(1,12),(1,15),(1,16),
          (11,13),(11,14),(2,15),(2,16),(10,13),(10,14)},
    H16={(0,10),(1,11),(1,14),(1,15),(10,12),(10,13),
          (2,14),(2,15),(9,12),(9,13)}.                      (29)

Every listed set is invariant under the corresponding labelled diagram's
automorphism group. Equivalently one may require every marked image
under every isomorphism with that standard target. This makes the
definition independent of an arbitrary identification of an actual
component with its labelled model.

### 8.5 The order12 and13 tables are short enough to print

At order13 the weights are all0 or1:

    A: unit exactly at cycle distance≤5;
    B: all pairs except(0,10),(10,11),(10,12);
    C: all pairs except(0,10),(0,11),(0,12),(1,12),(2,12);
    D: all pairs except {9,10}×{11,12}.                     (30)

Their costs are65,75,73,74. The D table has no half-weights at this
stage. In particular the later ambient58 center/fork refinements must
not be imported here.

At order12, before the ambient66-only patch, use

    A: unit exactly at cycle distance≤5;
    B: unit except weight1/2 at(9,10),(9,11);
    C: all pairs except(0,11);
    D: all pairs except {8,9}×{10,11}.                       (31)

Their costs are60,65,65,62. After the ambient66-only patch the costs
are60,64,63,62. Orders≤11 remain wholly paid throughout n66–97.

### 8.6 Finite data schema, exact reconstruction, and provenance

The mathematical supplement needs two finite kinds of table.

1. A preferred record is a tuple (L,M,f,a,b,J,p,q,h), where L,M
   are labelled minimal indefinite diagrams, f is an induced matching
   of overlapping nodes, a,b are their marked black pairs, J is the
   connected component at the selected corner, p,q are its two marks,
   and h≤17 is the proved maximum order of its whole completion.
   The routing appendix defines the exact admissible core domain and
   the deterministic selected-corner rule. Equivalent saved records
   contain the already extracted labels ofJ and the two positions.
2. A paid-set record is (b,m,X,P), where P is a list of distinct
   unordered pairs0≤i<j<m. It is interpreted literally as a unit
   indicator table. Half/zero patches are then(29)–(31), not inferred
   from costs. A finite rank17 record instead gives its labelled
   finite graph and a retained pair list T17(X).

For each preferred record with h≥m and |J|≤m enumerate every induced
label-preserving injection into the target, as described below. Unite
all its marked images with q_b. For P7 only order17 initially needs
this operation; its order16 successor uses all records with h≥16.
For P6 use m15–17; for P5 use m14–17. The applicable validated
removals are the nine cross-fork pairs X_m in D and the cycle pairs
above distanceb. Thus the precise finite operation is

    P_b(m,X)=(q_b ∪ all admissible preferred marked images)
               minus the stated validated removal set.     (32)

The injection algorithm assigns an unassigned pattern node and tries
every unused target node whose labels to all assigned nodes agree,
including zero labels. Necessary degree tests may reject a candidate,
but absence of an edge is never ignored. Its recursion depth is |J|≤17
and branching at most17, so it halts. Every induced injection has a
unique branch; induction on the assigned nodes proves completeness.
Deduplicating marked pairs takes their union, not their multiplicity.
Its mathematical correctness is separate from the geometric theorem
that these are the required corners and that each removed image has
a retained actual replacement.

For P5 the initial domain has949 overlap records,571 with h≥14.
The later distance-eight domain has2614 records,1153 with h≥14.
The latter's entire image union is checked to lie in the same paid
sets or their validated removal sets; it changes no table. For P6
there are437 eligible preferred records; for P7 order16 there are260,
and for order17 there are80. These numbers are consistency checks,
not the definition or completeness proof of the domains.

The rank17 finite table is reconstructed by taking the union of
R8's short/leaf9 pairs, the marked images of all bound17 selected
preferred records, and every marked image whose retained replacement
is a UNIT in the complete literal-interface ledger. Then remove A17's
endpoint pair(0,16), using its separate validated endpoint transfer.
It occurs in neither the basic nor the selected high-pattern union;
the complete literal ledger is recounted after this deletion rather
than assuming every earlier literal unit remains. This defines T17
without relying on a same-rank aggregate cost.
Its saved lists include every pair and complementary removed pair,
and the local proof validates all debited interfaces simultaneously.

For direct recovery of the finite supplement, the exact saved arrays
are the following. The common schema above, target graphs in8.1 and
construction(32) give their mathematical interpretation; no Python
function call is part of the rule definition.

|Table|Saved JSON array field|
|---|---|
|P7(17,X)|fresh_ideal7_finite_cross_payments: order17.targets[X].paid_pairs|
|P7(16,A/B/C)|fresh_order16_preferred_embeddings: targets[X].paid_pairs|
|P7(16,D)|fresh_order16_finite_cross_payments: exceptional_D15.remaining_D_pairs|
|P6(m,X)|fresh_ideal6_combined_payments: small_order_pair_sets, record keyed by(order,kind), paid_pairs|
|P5(m,X)|fresh_ideal5_small_factor_payments: records, keyed by(order,kind), paid_pairs|
|T17(X)|fresh_finite_d17_long_pair_support: records, keyed bykind, retained_pairs|
|Order13 C check|fresh_order13_payment_support: records[kind=C].required_union|

All stems in this table have suffix `_20261017`, and are in `results`.
Missing entries in an indicator paid-set mean zero. This convention
must not be confused with later `weighted_pairs` arrays, where one
checks their advertised complete pair domain explicitly.

The retained input file hashes, in table order, are

    ideal7_finite_cross: ae9b784fa99d7ee026bf792bab534d6810b5c633a13d910952088a053d1025ff
    order16_preferred:  28cbb8a2fcc0d13129f35e09e7c1a7bd6504d368aed130af33b70925e1738502
    ideal6_combined:    11bce448b82472ce582677e7d7609befa0886a8b19d54593905afb614c43e0bc
    ideal5_small:       be89d5f671965bccc75b9ea2ab8dff2a925b909b97a52e473c97a17bc85cdd34
    finite17_support:   a19e5b08125321f215b3a125225066376b6ee6dac3b91ad8f6b7a9d1fac8f328
    order13_support:    25b684c81c35e59c9b2afd2cb51f28d7847247b4b6d61e04c8670f8dbf494062

The D order16 retained list may equivalently be obtained by removing
X_16 from the preceding preferred table; this avoids a separate
mathematical table input. The source report contains the arbitrary-white
sealing proof for precisely that removal.

### 8.7 Local proof obligations and what remains

The local theorems establishing these rules are not just pair-count
lemmas. Their proof maps must show: every actual minimal core is in the
classified domain; every supporting-white contact is either excluded
by an exact Schur inequality, gives a retained short-path payment, or
seals the stated whole component; and every changed preferred or
literal interface still totals at least1. Triangle payments require
their separate retained-route theorem. No finite DP or cusp convolution
can substitute for any of these claims. The companion routing appendix
is responsible for these geometric derivations.

The actual rule evaluations for n66–97 are fully specified above.
To finish the full low tail, append the separate n59–65 patch ledger
starting from the ambient66 rule; its final ambient58 descendant is
not a valid substitute. The n59–62 exact varying-envelope values also
need to accompany that ledger. No additional computation is required
to read those finite arrays and overrides from their existing sources.

## 9. Source provenance and status

The definitions and exact outputs above were recovered by source reading
from these existing families: `fresh_ideal5_small_factor_payments`, the
`fresh_ideal5_weighted` and `fresh_ideal5_piecewise` interval certificates,
`fresh_ideal6_piecewise`, `fresh_order16_piecewise`, `fresh_ideal7_cross`,
`fresh_distance8_ideal`, `fresh_selective_distance9` through
`fresh_selective_distance13`, and `fresh_distance14_unrestricted`
(all20261017); and the ambient59–65 composition chain(20261018).
The geometry spine stages the exact interval endpoints. These names
are provenance, not mathematical definitions in place of Sections1–7.

The unbounded proof uses the existing
`fresh_nikulin_two_budget_frontier_20261015.py` and its independent
review, its dimension163 geometric binder, the independent
`fresh_nikulin_unit_charge_20261015_review.py`, and the grounded
`fresh_nikulin_height_charge_20261015_review.py` definitions. The
published narrow-facet theorem and finite/affine classifications are
external inputs. Framing, endpoint realization, common height, and the
compact-facet branch are explicit internal geometric inputs. Exact
KKT/packing checks establish the stated local numerical inequalities;
they are not a proof-assistant formalization or a replacement for those
geometric arguments.

This draft adds no theorem beyond the already retained proof chain and
makes no new execution claim. Its completion boundary is the separate
59–65 layer reconstruction and insertion of the full
local geometric lemmas supplied by the companion geometry/routing
appendix.
