# An A3 energy barrier and low-rank A-forest farthest vertices

20 October 2026. Paper proof only; no computation, enumeration, build,
optimizer or formal verification. Review is pending. Arithmeticity is
not assumed.

**Component theorem.** At a farthest vertex of a compact hyperbolic
Coxeter polytope whose whole finite diagram has only A1,A2,A3 components,
the pointwise energy in EACH A3 component is greater than1/4.

**Forest consequence.** In dimension n>=14, the farthest vertex cannot
have whole diagram A1^a+A2^r+A3^s. All other vertices are unrestricted.
This is a pointwise profile elimination, not a new universal compact
dimension ceiling.

The farthest-point and exact reciprocal-shadow mechanisms were developed
by root and batch_guard. They are rederived in the preceding source
`results/compact_a1_a2_farthest_vertex_energy_exclusion_20261020.md`,
SHA2563bcf38496cd0aa7ba922dd13101595a479c0c14937d5c08bb368e2df435707e2.
That source also proves the A1/A2 energy interfaces used below. The new
argument is the adjacent same-wall propagation in A3, its native
cost-gap obstruction, and the full audit of the additional A3-spectator
consumers. The pair/log-concavity shortcut and its1/4 application were
independently co-derived by root and resource_cleanup. No finite census
is imported.

## 1. Pointwise framework and propagation of the actual opposite wall

At an actual vertex p use p^2=-1, outward unit normals n_i, Gram B,
C=B^-1, d_i=C_ii and u_i=omega_i/sqrt(d_i). For a fixed interior future
unit point z choose p maximizing E=-z.p over the actual vertices. Set
x=-z_T/E, where z=E p+z_T. Then

    ||x||^2=1-E^-2<1,
    x.u_i>=tau_i=tanh(L_i/2).                           (1)

The energies ||x_component||^2 add over the orthogonal finite components.
For an outside actual wall c with height A=-p.c>0 and contacts
beta_i=-c.n_i>=0 put

    t=C beta,       v_i=t_i/sqrt(d_i),       T_i=v_i/A.

The exact actual-wall alternative is

    T_i>1: c itself closes edge i, tanh L_i=1/T_i;
    T_i<1: another wall closes it, tanh L_i>=T_i;
    T_i=1: impossible in a compact chamber.             (2)

For completeness, C-C_.i C_i./d_i is zero on row/column i and equals
B_(I\{i})^-1 on the other indices, hence is entrywise nonnegative.
For a DISTINCT opposite wall d with height D and contacts gamma, this
gives beta^T C gamma>=t_i(C gamma)_i/d_i. All-wall acuteness gives
A D>=beta^T C gamma, and the endpoint equation gives
tanh L_i=D sqrt(d_i)/(C gamma)_i. Thus tanh L_i>=T_i if c is not its
opposite. If T_i>1 that alternative is impossible, so c is its actual
opposite with reciprocal value. Value1 cannot give a finite endpoint
under either alternative. This proves (2) without an endpoint census.

There is a useful consequence for two directions i,j with positive
correlation k=u_i.u_j. If their actual opposite walls are DISTINCT,
apply the Schur inequality above to the opposite c of edge i. It gives
T_j>=k T_i=k/tanh L_i. Since c is not the opposite of edge j,

    tanh L_i tanh L_j>=k.                              (3)

For A3, in path order,

    C=[[3/2,1,1/2],[1,2,1],[1/2,1,3/2]],
    d=(3/2,2,3/2).

Adjacent direction correlations are k=1/sqrt3. Suppose this component
has energy at most1/4. With tau_i=tanh(L_i/2), equation(1) gives

    tau_i+tau_j<=sqrt[(1+k)/2]

for adjacent i,j. The function H(t)=2t/(1+t^2) is increasing on(0,1)
and log H is concave there: its second derivative is
-1/t^2-2(1-t^2)/(1+t^2)^2<0. Consequently

    tanh L_i tanh L_j
       <=H[(tau_i+tau_j)/2]^2
       <=32(1+k)/(9+k)^2<k.                            (4)

For the last strict inequality use k^2=1/3: after clearing positive
denominators it is equivalent to74k>39, and74^2>3*39^2.
Thus(3) forbids distinct actual opposite walls at adjacent indices.
Both adjacent pairs share their opposite, so ONE wall c closes all
three directions. Also(1) gives tanh L_i<=4/5 for every i. Therefore

    c closes all three directions,
    A^2<=(16/25)v_i^2 for each i.                      (5)

The remaining argument rules out every possible contact vector beta
in(5), retaining all spectator components.

## 2. Multiple contacts already exceed the small-energy window

Write the three contacts as u,v,w and let
q=beta^T C beta. The total squared height is A^2=q+Q-1, where Q>=0
is the cost of all spectator components. In particular A^2>=q-1.

Adjacent nonzero contacts are impossible: release the third old node,
whose actual opposite is c by(5). At that finite endpoint, the white
diagram contains c and the two adjacent contacted old nodes. Their
triangle is affine at the minimum contacts1/2 and is nonpositive if
either contact is larger, contradicting the positive-definite endpoint
Gram. This also excludes three contacts.

Otherwise a multiple contact must use only the two endpoints u,w.
If both are3, the old A3 path and c form the affine four3-edge cycle,
forbidden in a compact chamber. In every other case one contact is at
least1/sqrt2 and the other at least1/2. Then

    q-(16/25)v_2^2
       =(59/50)(u^2+w^2)+(9/25)uw
       >=177/200+9sqrt2/100>1.                         (6)

The last strict inequality follows from18sqrt2>23. Hence
A^2> (16/25)v_2^2, contradicting(5). This exhausts multiple contacts
with arbitrary additional spectator costs.

## 3. Exact spectator costs and the two single-contact cases

Because c closes an actual A3 edge, c plus the retained white normals
is positive definite. Hence on every WHOLE spectator component it has
positive-definite finite extension. For A1 the positive cost is a native
cosine square. For A2 and A3, two contacts in one spectator component
are impossible: adjacent simple3 contacts already give an affine
triangle; the two A3 endpoints give an affine four3-edge cycle.
Larger contacts only make these principal submatrices nonpositive.

The complete positive spectator costs below1 are therefore

    A1: cos^2(pi/m), m>=3;
    A2: 1/3,2/3,phi^2/3;
    A3: 3/8,3/4,3phi^2/8 at an endpoint, or1/2 at middle,
         where phi=(1+sqrt5)/2.                         (7)

For example the A3 inverse diagonal is3/2 at an endpoint and2 at middle,
so its allowed labels are respectively3,4,5 and3. Dotted or parallel
spectator contacts cannot have positive-definite extension.

The total Q is less than1. The following cost gaps follow directly from
(7), including sums over arbitrarily many components:

    Q>1/4  => Q>=1/3;
    Q>1/2  => Q>=7/12;
    Q>5/8  => Q>=cos^2(pi/5).                           (8)

Indeed three positive summands cost at least3/4. For two summands the
successive smallest sums are1/2,7/12,5/8,2/3; for one summand the listed
native gaps give exactly(8). Thus no unbounded finite-profile list is
hidden in these assertions.

Suppose c has just one ENDPOINT contact beta. By(5) c closes that
contacted direction; put r=tanh L_contact. Then

    r^2=1-(1-Q)/[(3/2)beta^2].                          (9)

For label3, Q>5/8 and(8) gives
r^2>=(sqrt5-2)/3>16/225. For label4, Q>1/4 gives r^2>=1/9.
For label5 the own-component cost is below1, so Q>=1/4 and
r^2>=sqrt5-2>16/225. Label6 is excluded by the affine3--6 path with the
next old node. For every larger finite label, and for a dotted contact,
beta^2>3/4 and (9) gives r^2>1/9. Thus ALWAYS

    r^2>16/225.

The inverse matrix gives T_far=T_contact/3. By(5) the same c closes
the far end, so tanh^2 L_far=9r^2>16/25, contradicting(5).

Suppose instead the only contact is at the MIDDLE. Then

    r^2=1-(1-Q)/(2beta^2).

For label3, Q>1/2 and(8) give r^2>=1/6. Label4 forms the affine
3/3/4 star with the two old endpoints. Label5 gives
r^2>=sqrt5-2>1/6; label6 gives an affine3--6 path; larger finite
or dotted labels give r^2>1/3. Hence r^2>=1/6 in every possible case.
Both endpoint shadows are T_end=T_middle/sqrt3 and are on the same
closed branch by(5). Their actual tangents therefore satisfy

    tanh^2 L_end=3r^2>=1/2.

The endpoint directions have inner product1/3. Their two constraints
alone require component energy at least

    (3/2)(sqrt2-1)^2=9/2-3sqrt2>1/4.                  (10)

The last inequality is289>288 after squaring17>12sqrt2.
This contradicts the assumed component energy. All contact vectors
have now been exhausted, proving the component theorem.

## 4. Audit of A3 spectators in the earlier A1/A2 energy bounds

Write a,r,s for the numbers of components. The preceding A1/A2 source
derives two interfaces, using thresholds
L*=asinh(1/sqrt2), sinh^2 L0=31/50 and e2=319/2000:

    E_A1>a/10-3/25;
    every A2 block pays>e2 except at most one isolated7 block,
    whose energy is>9/100;
    if isolated7 occurs, E_A1>=a/10.                   (11)

Here we check EVERY additional role introduced by spectator A3 blocks;
the old inventories and proofs remain as derived in the pinned source.

### A1-short directions

For a wall closing an A1 edge shorter than L*, with released contact
square b and white cost Q, the necessary open window is

    1-b<Q<1-(2/3)b.                                   (12)

Add the A3 costs in(7) to the old A1/A2 sums. For released label3 its
window(3/4,5/6) acquires no sum: new3/4 sums are boundary values, and
the next mixtures reach5/6 or exceed it. For released label4 the only
new sum is A3_end3+A1_3=5/8; it contains the affine3/3/4 star.
For released label5 the new costs are A3_end3=3/8 and A3_middle3=1/2.
For released label6 only A3_end3=3/8 is new, and it gives the affine
3--6 path. For finite labels>=7 the upper bound is below1/2, so the
only new cost is a SINGLE A3_end3. It is allowed by(12) only if
cos^2(pi/m)<15/16, hence m<=12; label13 already exceeds15/16 since
sin(pi/13)<22/91<1/4. Dotted releases have upper bound below1/3.

Thus all additional A1-short walls are precisely the POSSIBLE profiles

    P=(A1_5,A3_end3),       height squared (sqrt5-2)/8<1/32;
    Q5=(A1_5,A3_middle3),   height squared (sqrt5-1)/8;
    R_m=(A1_m,A3_end3),     m=7,...,12,
                           1/8<height squared<5/16.   (13)

Each such wall has only ONE A1 direction and cannot coexist with any
other A1-short wall. Here are the full compatibility checks.
All A1-short wall heights are below1/2, so distinct ones are orthogonal.
Every positive overlap of their spatial supports costs at least1/8,
because the smallest inverse entry among A1,A2,A3 is the A3 endpoint
entry1/2 and each contact is at least1/2. Therefore P cannot coexist
with any other: its height product is below1/8, and disjoint support
would give a forbidden product strictly between-1/2 and0.

For Q5, a shared A1 contact contributes at least phi/4, too much. A
shared A3 endpoint contributes1/4 while the other new wall's squared
height is below5/16; their height product is less than1/4. Sharing its
A3 middle contributes1/2. Old A1-only or A1/A2-mixed short walls have
no A3 overlap. These alternatives exclude Q5 with every other profile.

For two R walls, a shared A1 contact is again too large. In the same
A3 block, equal endpoints contribute3/8, greater than their height
product; opposite endpoints contribute1/8, LESS than their height
product (each squared height exceeds1/8). Neither is orthogonal.
Different components give no overlap. With an old A1-only or
A1/A2-mixed short wall there is either no overlap, or an A1 overlap
at least sqrt3/4, greater than their height product (<sqrt(5/32)).
The P,Q5 cases were already excluded. This proves the assertion.

Consequently the A1 debit in(11) is unchanged. If a new profile occurs
there is at most one short A1 direction, whose loss from the1/10
baseline is less than1/10<3/25. Otherwise the complete old argument
applies, including its two-direction and energy estimates.

### A2-short directions

The old proof classifies a contacted A2 direction shorter than L0 by
its released square b and off-block cost Q:

    max(0,1-(4/3)b)<Q<1-(200/243)b.                    (14)

For released label3, the additional A3 sums in its window
(2/3,193/243) are

    3/4: one A3_end4, two A3_end3, or A3_middle3+A1_3;
    17/24: A3_end3+A2_3.

An A3_end4 gives a principal affine chain3--3--4--3. The two-A3 and
A3/A2 cases give a principal affineE6 with three arms of length2.
The A3_middle3+A1_3 case gives the affine double-fork diagram with
two adjacent degree3 nodes. For label4, the new sums in
(1/3,143/243) are3/8 and1/2. The first contains the same affine
3--3--4--3 chain. For the second, take c, the A2 contacted root,
and the A3 middle with its two leaves: these have positive null
coefficients(1,1/sqrt2,1,1/2,1/2), respectively. Thus it is also affine.

The named affine checks require no table: the3--3--4--3 chain has
positive unit-Gram null coefficients(1,2,3,2sqrt2,sqrt2); the E6 star
has central coefficient3 and coefficients(2,1) on each arm; the
double fork has coefficients2 on both branching nodes and1 on its
four leaves. Direct substitution gives zero in each case.

For label5 the only new sum is3/8, giving

    G=(A2_5,A3_end3),
    H_G=height squared=phi^2/3-5/8<1/4.                (15)

Label6 remains affine-forbidden, and labels>=7 have upper bound below
1/3, less than the least A3 cost3/8. These exhaust all additional sums
in(14): another positive contact costs at least1/4, and the ordered
costs are exactly(7).

Profile G cannot close its uncontacted A2 mate, since that white cost
is cos^2(pi/5)+3/8>1. Apply the reciprocal shadow to its ACTUAL distinct
opposite. It gives

    tanh^2 L_mate>=cos^2(pi/5)/(3H_G)>2/3>16/25,

so this block has energy>1/4>e2. The A2 accounting in(11) is unchanged.

Finally an isolated7 A2 wall remains incompatible with EVERY A1-short
wall: the old profiles are excluded in the pinned source, while each
new profile in(13) has support only in A1/A3 and hence disjoint from
that A2 wall. Their height product is below1/2, in the forbidden gap.
Therefore the nonconcurrent A1 and isolated7 debits in(11) also remain
valid. No inherited coupon is silently discarded.

## 5. The full mixed forest has more than one unit of energy

Let n=a+2r+3s>=14 and set e3=1/4, e2=319/2000. The component
theorem and(11) give, if isolated7 occurs,

    ||x||^2>a/10+e2*r+e3*s-7/100.

Since e3/3>e2/2 and1/10>e2/2, the right side is at least

    (319/4000)n-7/100
       >=(319/4000)*14-7/100=2093/2000>1.             (16)

If isolated7 does not occur and a>=1, the bound is

    ||x||^2>a/10+e2*r+e3*s-3/25
             >=(319/4000)n+(81/4000)a-3/25
             >=4067/4000>1.                           (17)

If a=0 there is no A1 debit at all. The A2 and A3 blocks instead
contribute at least(319/4000)n>=2233/2000>1. Strictness is supplied
by every present component; absent component counts incur no fictitious
strict credit.

All cases contradict ||x||^2<1 in(1). This proves the forest consequence.

The proof has not controlled A4 or larger components, or B,D,E,F,H
components. Their possible occurrence at the farthest vertex remains
the obstruction to a global dimension ceiling. The new component gap
and its branch-propagation mechanism supply a concrete uniform route
past isolated continuation counts, rather than a classification of all
compact polytopes or all arithmetic fields.
