# Mathematical algorithm fragment for the finite-volume dimension bound

20 October 2026. Staging text for the existing `finite_volume57.tex`;
this is not another paper. No producer, replay, search, Lean or TeX job was
run in preparing it. The primary manuscript and frozen evidence were not
edited. The results below describe the existing computations and give
their mathematical correctness arguments. They do not supply the separate
geometric routing and whole-component theorems.

## A. Finite forests: the mathematical input and bounded domain

The finite calculation concerns rank58. Multiply every charge by64,
so all costs in this subsection are integers. The geometric local rule
supplies two finite tables: the ordinary cost `b(t)` of a connected
finite type, and the symmetric ordinary interaction `j(t,u)` between
two distinct components. Thus `j(t,t)` concerns two copies of the same
type; it does not repeat a charge within one component. Presence-dependent
costs are kept separately in the function `Psi` below.

Use the finite alphabet

    T = {A_r:1<=r<=58} union {B_r:2<=r<=58}
        union {D_r:4<=r<=58}
        union {E6,E7,E8,F4,H3,H4,I2*}.

Here `rank(I2*)=2`. The last symbol represents every remaining finite
dihedral type I2(m), m>=5, for this cost calculation. It is justified by
the all-label local-rule theorem, which proves that these types have the
same relevant charge and interactions. It is not a truncation of m and
does not identify their metrics. The cost of a forest with multiplicities
`m_t` is

    C64(m) = sum_t b(t)m_t
             + sum_t j(t,t) binom(m_t,2)
             + sum_(t<u) j(t,u)m_t*m_u + Psi(m).             (A1)

The domain is exactly

    m_t integer >=0,   sum_t rank(t)*m_t=58.                 (A2)

Every count is at most58/rank(t), so (A2) is finite. An actual finite
vertex has a forest in this domain by finite Coxeter classification.
The converse is unnecessary: formal forests need not be realizable.

For clarity the complete context function of the *preceding* rule is
given here explicitly. Write `a_i=m_(A_i)` and `d_i=m_(D_i)`, and let
`[condition]` denote0 or1. In the same integer scale,

    Psi = 36*a2*(a2-1) [a1>=2 and a10>0]
        + 64*a2*d5 [a1>0 and a9>0]
        + 24*a1*d5 [a11>0 or (a9>0 and a2>0)]
        + 40*a1*d5 [a9>0 and a2>0]
        + 32*(a15+d15) [a1>0]
        + 64*d14 [a1+a2+a3>0]
        + [a1>0]*(80*a16+32*a14+64*d16+64*a13
                     +32*a4*a12+64*a2*a4*[a10>0]).          (A3)

In particular witness presence is Boolean: a second A1 or A10 does not
multiply a presence bonus, except where the displayed formula explicitly
contains its multiplicity. Formula (A3) is the flattened mathematical
version of the contextual costs used by the final preceding finite
count. The ordinary tables exclude these terms.

**Data/proof boundary.** To make the paper self-contained, its supplement
must print or mathematically define the complete tables b,j and prove
that (A1) is the charge of the simultaneous local rule. A program import
chain is not a substitute for these two definitions. The recurrence
below is fully specified for any supplied b,j,Psi, and its application
uses the frozen tables; the geometric justification of those tables is
a separate lemma. This fragment deliberately does not claim to have
re-proved that universal lemma by reading table entries.

## B. A residual maximization and its independently checkable certificate

Fix a finite residual alphabet U, positive integer ranks r(t), integer
weights w(t), and a remaining rank R with0<=R<=58. Define

    V(0)=0,
    V(s)=max { V(s-r(t))+w(t) : t in U, r(t)<=s,
                                  V(s-r(t)) is reachable } (B1)

for1<=s<=R. An empty maximum is marked *unreachable*, not assigned an
ordinary integer value. Evaluate (B1) in increasing s. Save one last
item attaining the maximum at each reachable positive s.

**Lemma B1 (termination and exactness).** This algorithm halts. Its value
V(s) is exactly the largest total weight of a multiset of U of rank s;
an unreachable state means that there is no such multiset.

*Proof.* Every recursive predecessor has smaller rank, since r(t)>0.
There are at most59 states and finitely many transitions per state.
For induction on s, removing any item from a nonempty multiset gives
one of the candidates in (B1), proving its weight is at most V(s).
Conversely an attaining item appended to the induction's attaining
multiset realizes each finite candidate. If no candidate exists, a
nonempty multiset would supply one by removing an item, a contradiction.
The empty multiset gives the base case. Repeatedly removing a saved
last item terminates after at most58 removals and recovers an attainer.

The saved-table checker need not rerun the maximizing recurrence. It
checks V(0)=0 and, at every positive state s:

1. For every reachable predecessor, V(s) is reachable and
   `V(s)>=V(s-r(t))+w(t)`.
2. If any reachable predecessor exists, the saved last item is such a
   predecessor and its inequality is an equality.
3. If none exists, V(s) and the last item are both unreachable.

The same induction proves exactness from these three checks. The upper
inequalities alone prove an upper bound but not attainment or correct
unreachability. Equality witnesses without *all* upper inequalities do
not prove that an omitted multiset has lower cost.

An implementation may generate the same values by an unbounded
item-first recurrence. Its checker above uses every possible last item,
so correctness does not depend on trusting that generation order.
Machine sentinels are only encodings: the mathematical value of an
unreachable state remains separate from every integer.

## C. Context reduction used for the finite rank58 calculation

The calculation first bounds the domain by representatives, and then
refines every class that can exceed377. Both steps are needed; refining
only previously recorded maximizing forests would not be complete.

The representative map sends every rank6 type to E6, every rank9 type
to B9, A12 and B12 to A12 while keeping D12 distinct, and keeps A17
separate. All remaining ranks at most8 map to A of that rank; all other
ranks map to D of that rank. In particular B17 and D17 remain residual
D17 items rather than joining the fixed A17 count.

The exact outer domain is

    c9,...,c17 >=0,  sum_(r=9)^17 r*c_r<=58,
    0<=a12<=c12,
    0<=k5<=floor((58-sum r*c_r)/5).                          (C1)

All variables are integers. The fixed representative multiset has
`c9` B9, `c10` D10, `c11` D11, `a12` A12,
`c12-a12` D12, `c13` D13, `c14` D14, `c15` D15,
`c16` D16, `c17` A17, and `k5` A5 components. Its residual ranks are

    {1,2,3,4,6,7,8} union {17,...,58}.                      (C2)

At fixed components F, a residual actual type t has ordinary item weight

    w_F(t)=b(t)+sum_(f in F) j(t,f).                        (C3)

All ordinary interactions between the allowed residual types are zero;
this is checked for the complete finite alphabet, not inferred from an
attaining forest. Context-sensitive counts are fixed or separated into
the branches below before (B1) is applied.

The inherited representative cover has1295 fixed large vectors after
the A/B12 versus D12 split and3298 expanded choices including k5. For
each expanded choice it has an ordinary and a no-A2 residual table,
6596 in total. The interaction key is

    (2c10+3(c11+...+c16)+2c17, c11+c12, c13, c11,
                       c14, c9, c10, c16, c15, a12).       (C4)

This key retains all information used by those representative item
weights. Rank5 is already fixed; it is not present in (C2). The
representative-domination checks bound every actual base and cross
interaction, retaining the explicit22 per fixed A17 increment. The
context terms not dominated by this cover are exactly the two families
described below. Outside them, all3298 expanded rows are checked, and
exactly39 have upper value greater than377.

### C.1 All actual fibers and exact residual contexts

For each of these39 classes, expand every fixed representative component
into its full actual-type fiber. Identical components are expanded by
nondecreasing type choices, so there is no missing mixture or unnecessary
ordering distinction. This gives748 fixed actual multisets F. For one
such F set R=58-rank(F) and explicitly enumerate

    0<=a1<=R,  0<=a2<=floor((R-a1)/2).                      (C5)

These are the complete actual A1 and A2 counts. Neither type remains in
the residual alphabet. If a1,a2 and a fixed A10 are all present, exclude
A4 from that alphabet: all such forests containing A4 belong to the
separately bounded family in C.2. Otherwise A4 is allowed.

The residual context flag that is not already fixed by (C5) is the D14
bonus triggered by an A3 when a1=a2=0. If F contains D14, split this
case into:

* no A3: remove A3 from the residual alphabet;
* at least one A3: reserve one A3, charge its ordinary weight and the
  single bonus64 times the number of D14 components, and allow further
  A3 items in the residual table without repeating this presence bonus.

In all other cases no reservation is needed. With `has1=[a1>0]`, a
residual A4 item also receives `32*has1*m_F(A12)` from (A3). All other
terms of Psi are determined by F,a1,a2 and the reserved-A3 branch.

Write eps=0 or1 for the reserved A3 and let U be the permitted residual
alphabet. The candidate is

    cost(F)+a1*w_F(A1)+a2*w_F(A2)+Psi(F+A1^a1+A2^a2)
       + eps*(w_F(A3)+64*m_F(D14))
       + V_U(R-a1-2a2-3eps),                              (C6)

where `cost(F)` in (C6) means its *ordinary* cost. The weights in V_U
include the just-specified guarded A4 increment. An unreachable residual
rank contributes no candidate. The two A3 modes are used only when
a1=a2=0 and D14 is present; in the reserved mode the displayed Psi
therefore has not already charged the A3 bonus. Thus no presence charge
is omitted or paid twice.

Take the maximum over (C5) and its allowed modes. Each forest in the
complement determines one candidate of (C6): remove its fixed types,
then its A1/A2 factors, and reserve one A3 exactly when required.
Conversely reconstruction from a reachable residual state gives a
forest in the complement with the claimed full contextual cost. This
proves both completeness and exactness on every actual fiber.

All residual table keys are `(has1, excludeA4, excludeA3)` at fixed F;
there are at most8 keys per fiber and at most59 rank states per key.
Every outer count and every type-fiber choice is bounded by (C1),(C5).
Hence the full context algorithm halts, even without using the smaller
recorded counts. Certificate replay checks the table inequalities of
Lemma B1, then reconstructs each winning forest and evaluates (A1)
directly. It does not infer contextual correctness from additive
attainment alone.

### C.2 The two separately bounded families

The first family is the extra guarded A1/D5 charge under A9 and A2.
The inherited complete759 conditional rows include this charge. Their
representative allowances, with22 added for every fixed A17, bound this
family by2837/8<377. This is an upper relaxation, not an attainment
statement.

The second family has a1,a2,a4 all positive and a whole A10. At each
fixed representative/rank5 choice of remaining rank R, enumerate

    1<=a1<=R-6,
    1<=a4<=floor((R-a1-2)/4),
    1<=a2<=floor((R-a1-4a4)/2).                            (C7)

Set L=R-a1-4a4-2a2. Use the inherited no-A2 table at rank L, adding
the fixed ordinary costs and all contexts evaluated from the explicit
counts, in particular64*a2*a4 and the A17 increment. Extra A1 or A4
allowed by this upper table only enlarge its domain. This is therefore
a valid coarse upper bound, with539 classes and37,985 count choices.

There are14 classes above377 in this coarse relaxation. Refine *all*
of them into the77 actual fixed assignments containing A10, enumerate
(C7) again, and use their actual residual tables. Their weights dominate
all remaining ordinary items and fixed interactions; no remaining item
can activate an omitted context. The residual relaxation may still
allow extra A4. All23,376 choices have value at most2861/8<377.
These inequalities also cover overlap with the first increasing family.

### C.3 Computed finite proposition and the final monotone change

The following is the mathematical result recorded by the existing
complete count and separate saved-table replay; it is not a new run:

**Computed proposition C.** For the preceding simultaneous rule, every
rank58 forest has charge at most377, except exactly the eight forests

    B13 or D13, together with three factors from {A15,D15}.

Each of the eight has charge378. The complementary calculation has748
actual fibers,1549 residual tables,96,528 contextual choices,758,838
last-item upper inequalities and31,545 equality predecessors. In the
only class with value above377 the residual rank is zero, so its eight
listed actual fibers exhaust the obstructions, rather than merely
exhibiting eight maximizing witnesses. The cost there is78+3*100.

The final geometric rule decreases B13 by1. Its additional contextual
change, in the same integer scale, is

    -64*d13*[a15+d15>=3]
       + (32*d28+64*a12*d16)*[a15+d15>=2].                 (C8)

The guards concern distinct whole components, as proved geometrically.
Rank58 capacity makes every increasing forest exactly D28 plus two
rank15 factors, or A12+D16 plus two rank15 factors. There are three
multisets of two A15/D15 factors, hence six increasing forests. Their
new charges are respectively370.5 and375.5. The eight previous
obstructions drop to377, and outside these six families the change is
pointwise nonincreasing. Thus the final exact *formal* maximum is377.
No new Bellman calculation is needed for this last step, and no
geometric realization of a formal attainer is claimed.

## D. Terminal cusp enumeration: an independent finite algorithm

The link dimension is d=57. A cusp's affine factor ranks are positive
integers r_i summing to57. The local factor-fiber theorem supplies the
normalized envelope c(r). Its complete numerical input can be printed
compactly as follows:

    c(r)=r/2                           (1<=r<=8),
    c(r)=(5r+56)/(r+1)                 (24<=r<=57),

and

|r|9|10|11|12|13|14|15|16|
|---|---|---|---|---|---|---|---|---|
|c(r)|71/16|53/11|31/6|11/2|40/7|367/60|411/64|227/34|

|r|17|18|19|20|21|22|23|
|---|---|---|---|---|---|---|---|
|c(r)|119/18|128/19|137/20|146/21|155/22|163/23|85/12|

These are all required ranks, not a sample. The all-family geometric
envelope lemma proves `C_actual<=V*sum c(r_i)`; the arithmetic algorithm
does not establish that lemma.

Let H(n,l) be the list of nondecreasing partitions of n with parts at
least l. Define it recursively by

    H(0,l)={()},
    H(n,l)=union_(q=l)^n {(q) concatenated with p:
                                      p in H(n-q,q)}       (D1)

for n>0, with an empty list when l>n. The terminal domain is exactly

    T57 = union_(ell=0)^7 { (ell) concatenated with p:
                              p in H(57-ell,8) },         (D2)

where the zero part is omitted. Thus every part is at least8, with at
most one positive exceptional part1,...,7.

**Lemma D1 (finite complete enumeration).** Recursion (D1) halts and
(D2) lists each terminal partition once.

*Proof.* Along a recursive call the first argument decreases by the
positive integer q. Every loop has at most57 choices. A nonempty
nondecreasing partition has a unique first part q and a unique tail
with parts at least q; induction proves both coverage and uniqueness
in (D1). The possible exceptional part below8 is unique, so the eight
branches in (D2) are disjoint. A terminal partition without such a part
belongs only to ell0. This proves the assertion.

One may independently enumerate multiplicities k8,...,k57 with
`0<=kq<=floor(remaining/q)`, accept remaining rank0,...,7 as the
exceptional part, and compare the two exact lists. This uses a different
branch decomposition; it is a useful implementation check, not an
additional mathematical hypothesis.

For a rank list r, multiply the factor polynomials by exact convolution:

    H0(t)=1,
    H_(k+1)(t)=H_k(t)*(1+t+...+t^(r_(k+1))).                (D3)

Equivalently set `new_h[j]=sum_(q=0)^r h[j-q]`, with out-of-range
coefficients zero. After all factors, the degree is57. Compute

    V=prod_i(r_i+1),
    G=sum_(j=0)^57 min(j,57-j)*h_j,
    kappa=(5G-57V)/2,
    m=kappa/V-sum_i c(r_i),  B=V*m.                        (D4)

All arithmetic is integer or rational. Every polynomial has58
coefficients. The recursion depth is at most57 (at most7 high parts
in the terminal domain), and (D3) is a finite nested sum. These explicit
bounds prove termination of the entire arithmetic calculation.

Checking only polynomial degree, symmetry and coefficient sum would not
verify G. Recompute every coefficient by (D3), then compute (D4).
Compare the output list of rank tuples with the independently generated
domain, rather than trusting a saved number of rows.

**Computed proposition D.** The existing complete enumeration and
independent endpoint replay give exactly2552 terminal profiles, and

    min_(r in T57) m(r)=87/10240>0,
    min_(r in T57) B(r)=36.                                (D5)

The two minima need not have the same attaining profile; the proof
uses each for its stated purpose. These are assertions about the
complete domain (D2), not sampled cusp shapes.

## E. Why terminal verification bounds every cusp

For1<=a<=b<=7 replace the two ranks(a,b) by(a-1,b+1), omitting zero.
This preserves total rank, and uses only c(r)=r/2 through rank8, never
the changed rank9 entry. It therefore preserves the normalized cost.
The integer potential sum_i r_i^2 increases by

    (a-1)^2+(b+1)^2-a^2-b^2=2(b-a+1)>=2.                  (E1)

It is bounded above by57^2 and initially at least57. Thus at most1596
transfers occur. If two positive ranks below8 remain, another transfer
is available. The process therefore halts in exactly the domain(D2).

Here is the finite certificate for the curvature comparison. Let
U_a,U_b be independent uniforms on{0,...,a},{0,...,b}. Define

    P_(a,b)(j)=#{(u,v):u+v=j}/((a+1)(b+1)),
    Delta_t=sum_j max(j-t,0)*(P_(a-1,b+1)(j)-P_(a,b)(j)).  (E2)

Check total mass equality, equal means and Delta_t>=0 for
t=0,...,a+b, for all28 pairs1<=a<=b<=7. A stop-loss difference
is linear between consecutive integer thresholds and zero beyond the
finite supports, so these checks give convex order. Equivalently,
every convex function on this finite interval is an affine function
plus a nonnegative combination of its integer stop-loss functions.
The existing certificate verifies all28 exact rational lists.

Condition on all other independent uniform factors to see that adding
them preserves convex order. In particular the transfer increases
`E|2S-57|`. Since

    G/V=(57-E|2S-57|)/2,

it does not increase kappa/V and hence does not increase m. At the
same time the vertex product decreases, because

    (a+1)(b+1)-a(b+2)=b-a+1>0.                            (E3)

Writing original and terminal quantities as(V,m) and(V',m'), one has
`V>=V'` and `m>=m'>0`. Consequently

    kappa-C_actual >= V*m >= V'*m' >=36.                  (E4)

The strict positivity from (D5) is essential: multiplying the two
separate inequalities would not be valid for a negative terminal
margin. This proves the claimed ideal surplus for every actual cusp,
including those with several affine factors, once the factor-fiber
geometric envelope is established.

## F. Translation and execution provenance

The mathematical algorithms above have these existing implementations:

|Statement or layer|Source location|
|---|---|
|Representative domain and two table modes|`fresh_a16_d16_witness_context_ambient58_frontier_20261018.py`, `finite`, `table`, `check_table`|
|Rank12 type split and fixed fibers|`fresh_d12_type_split_ambient58_forest_frontier_20261018.py`, `representative`, `fixed_domain`, `fixed_components`|
|Final preceding context polynomial|`fresh_d17_leaf_a16_context_coupled_transfer_20261018.py`, `actual_base`, `actual_cross`, `context`, with its pinned predecessor definitions|
|Actual complement, guard modes and last-item replay|`fresh_d17_a16_coupled_ambient58_frontier_20261018.py:105-185`|
|Both exceptional-family bounds|Same source, `representative_cover` and `increasing_family`|
|Complete low8 domain, independent reconstruction and envelopes|`fresh_order10_low8_ambient58_frontier_20261018.py`, `high8_parts`, `concentration`, `ideal`|
|Independent2552 convolutions and28 stop-loss replays|`fresh_d13_context_dimension57_20261018.py:33-92`|
|Final eight-plus-six monotone finite closure|Same endpoint and `fresh_d13_three_rank15_context_transfer_20261018.py`, `finite`|

The frozen complete preceding certificate and its report record successful
generation and separate saved-table replay of Proposition C. The frozen
low8 certificate records generation and replay of Proposition D, and
the final dimension57 endpoint independently re-enumerates(D2) and
recomputes(D3)--(D5) and all28 transfers. The endpoint does **not**
rerun the earlier finite Bellman domain or prove the geometric atlas.
No process identifiers or historical debugging transcript are part of
these mathematical statements.

## G. What this does not yet make standalone in the paper

These are missing *expositions*, not a declaration that the internal
theorems have no proofs. Their source proofs remain necessary.

1. **Ordinary weight data and universal routing.** Print the complete
   intrinsic weight/base/cross tables and prove their simultaneous
   validity for arbitrary supporting whites and arbitrary native labels.
   This includes actual-face support, the high-core selected/nonselected
   split, whole-component seals and chosen-output retention. A list of
   saved small cores cannot replace the no-missing-support theorem.
2. **Why all finite cases enter the context cover.** Give the mathematical
   representative-domination inequalities and the exhaustive separation
   of the two increasing families. The formulas above specify the
   algorithm, but do not prove that a geometrically valid payment is
   dominated by its tables. The presence of arbitrary I2(m) is a
   quantified local lemma, not a finite-label assumption.
3. **The dimension58 final exchange.** The written selected33-context
   split, arbitrary indirect-white rank bound, affine-D12 quiet and
   energy-one repairs, B13 endpoint consumers and retention of the
   actual triangle payments must be connected into one simultaneous
   rule. The current manuscript contains much of this endpoint proof;
   its older routing inputs still reside in the atlas.
4. **Ideal geometry.** Prove the full nonsimple cusp factor-fiber
   multiplicity and whole-affine-factor containment, then derive the
   numerical envelope c(r) from the same simultaneous rule. Polynomial
   enumeration cannot establish those geometric facts.
5. **Global shape and higher tail.** State the exact published
   incidence/index identities with hypotheses, prove the internal
   two-extrema correction, and expose the completed bound58 argument
   excluding all dimensions at least59. Excluding dimension58 alone
   gives no higher tail, and general Coxeter facets cannot be assumed
   to supply it by descent.

The companion lemma map assigns these obligations explicitly. Integrating
this algorithm fragment closes the definition/termination/completeness
gap for the arithmetic steps. It does not convert the current concise
endpoint paper into a self-contained proof of the entire geometric
atlas or its higher-dimensional tail.
