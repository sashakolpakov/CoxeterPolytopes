# A19 and A20: actual projection above six sevenths

20 October 2026. Paper reasoning by root, held for full outside-author
review. This is a supporting research note, not a new manuscript. No
scientific computation, optimization, Lean, TeX or accepted-source edit.

**Theorem.** At a vertex farthest from a strict interior point of a
compact hyperbolic Coxeter polytope, an A19 or A20 component has

    ||x_A||^2 > 6/7,

with arbitrary finite spectators. The conclusion is about the ACTUAL
projection of the common farthest vector, not an edge-only relaxation.
No new global ceiling is asserted before the other dimension23 cases
are separately proved and the complete forest accounting is reviewed.

The argument needs neither the weak-index-five exclusion nor the long
plateau-energy certificate used to exclude A21. At this lower target,
two separated strong directions already suffice. We do need the complete
weak-profile inventory, including its actual-point and spectator clauses.
The rank extension of those clauses is checked explicitly below.

## 1. Exact internal inputs and the smaller rank range

All paths in this note are under results/.

| Input | Source and SHA256 |
| --- | --- |
| Complete central finite-white inventory and isolated-cap estimates | `compact_whole_a21_farthest_vertex_exclusion_candidate_20261020.md`, `3abe11bb70ce677069279462a91111ed6b8ab8889e3ec9e656ac71c1094200fb` |
| General-N adjacent-packet algebra and exact mixed-pressure bounds | `compact_whole_an_geometric_a21_plateau_transfer_20261020.md`, `8a58b7747e1a3d8c9affe7df17038781a62983ddd70ccb670502c6f40b1f9cee` |
| Complete zero-own-contact repair and clearing arbitrary spectator pressure | `compact_an_large_components_arbitrary_spectator_transfer_20261020.md`, `8d3b999761be7dd46e7990973b3ec4d4a746d4d682b044d53695f4aae405e3e5` |
| Actual mixed finite replacement, including native products and full displacements | `compact_a21_adjacent_weak_far_m_simple3_packet_20261020.md`, `1c4437c98d916d558081572dd2ecc12bc75eda1328a53633fe111120ec63dcb7` |
| Far5 distance-one prohibition and actual near4 mate, valid N>=15 | `compact_an_projection_pressure_and_weak_cap_packing_review_20261020.md`, `8b0a59493f1fa5dc7081e5e70beff2839ac5918bfcd50e089e05ae9b5546e43b` |

The published inputs are finite spherical classification, compact acute
simplicity and the native gap for actual walls. The listed new geometric
lemmas are internal proved inputs, not published assumptions. We use the
proofs of their inventories and pressure constructions, rather than
importing their rank21 conclusions at an unproved smaller rank.

Put N=r+1 in {20,21}, c0=61/100, and use unit A_r normals
n_i=(e_i-e_(i+1))/sqrt2. Their spatial duals have

    D_i=omega_i^2=2i(N-i)/N,
    C_ij=omega_i.omega_j=2i(N-j)/N, i<j.

Let h_i=x_A.omega_i, h_0=h_N=0. The strict interior point is
E(p-x), and

    e=||x_A||^2<1,
    e=(1/2)sum_(i=0..N-1)(h_(i+1)-h_i)^2.

The h_i are strictly concave because x_A.n_i>0. On the central
window 5<=i<=N-5, call i weak if h_i<=q_i-c0, q_i=sqrt(D_i).
The exact bounds needed in the extended inventory are

    D_i>=15/2, q_i>27/10, q_i<31/8,
    C_(i,i+1)>=7.

We argue under the contrary assumption e<=6/7 throughout.

## 2. Complete weak inventory, including arbitrary spectators

At the actual opposite c of a central direction i, the retained old
A paths have ranks i-1 and r-i, both at least4. Finite spherical
classification therefore gives exactly the same four alternatives as
in the pinned whole-A proof:

* no old A-path contact;
* one simple3 contact at an admissible end, D position or E position;
* one end4 contact;
* two end3 contacts, one to each path.

H/F completions are too small. In the last alternative the connected
white completion has rank r>=19>8; the arm and bond arguments in the
pinned source force A_r. No step here requires r>=21.

Write c=A p-v_A-v_S, Q=||v_S||^2. The actual endpoint retains the
entire spectator forest, so Q<1. A contacted spectator contributes
at least1/4. The zero-own-contact repair in the pinned component
transfer uses only the retained path length m>=5 and old affine
principals, which exist already at r>=19. Its remaining far-end3,
next-to-end3 and far-end4 branches use respectively Q>=1/4,
Q<2/m and Q<1/m to give either negative height or an actual
orthogonal opposite. These statements and their contradictions are
unchanged. Thus every central opposite has positive own contact.

For an isolated own contact gamma, positive clearing gives the exact
same pure-component bounds as in the pinned transfer. With
u=(1-Q)/gamma^2, they are

    h_i > sqrt(D_i-1)-sqrt(u-1),                 if u>1,
    h_i > sqrt(D_i-u)-sqrt(u(1-u)),             if 0<u<1.

The case u=1 gives a feasible nonzero null ray and is impossible.
No cleared normal is assumed to be an actual facet.

Here are the only scalar changes needed at N20/21. From D_i>=15/2,

    q_i+sqrt(D_i-1)>27/10+5/2=26/5.

Consequently q_i-sqrt(D_i-1)<5/26. For Q>0 and u>1 the same
native bound gives sqrt(u-1)<5/12, so the total loss is

    q_i-h_i < 5/26+5/12=95/156<61/100.

For 0<u<1, either with spectators or for a dotted isolated cap,

    q_i-h_i < (5/26)u+sqrt(u(1-u))
              <=(5+sqrt701)/52<63/104<61/100,

using sqrt701<53/2. For Q=0 native m>=8, the finite-face bound
uses tan(pi/m)<=tan(pi/8)<5/12 and the same95/156 seal.
Native7 uses q_i>5/2 and sqrt(D_i-100/81)>5/2; the latter follows
from 15/2-100/81>25/4, with difference5/324. Its exact loss is
therefore less than

    20/81+9/25=1229/2025<61/100.

The isolated5 proof is unchanged: q_i+sqrt(D_i-25/16)>19/4
and its loss is less than25/76+5/24<61/100. Labels3,4,6 and
parallel contacts retain the same affine exclusions. Thus every
isolated cap is strong at these smaller ranks as well.

For the other finite-white cases, the edge-target ratio is still

    F(q)/q < 761/500,
    F(q)=[1+(1-c0/q)^2]/[1-(1-c0/q)^2],
    5/2<=q<=31/8.

This is the positive concave-quadratic seal in the pinned general-N
proof. All nonexceptional role bounds use only m>=5: the near3
native>=5 bound (101/125+2/5)/sqrt(3/5), the D-position bound
(101/125+1/5)/sqrt(2/5), and the larger E/end4/two-white3 bounds.
They exceed761/500. The same local affine exclusions are principal
subsets of the actual full star. Adding foreign cost Q strengthens
these edge inequalities. The only exceptional near4/simple3 and
far-end3 roles with Q>0 are made strong by the exact estimates in
Section4 of the pinned component transfer; those estimates also
use only m>=5 and the displayed F(q)/q bound.

We have therefore obtained the COMPLETE weak list:

* the literal near4/simple3 cap, with Q=0;
* the literal far-helper profile

      v_A=gamma*omega_i+omega_e/2,
      e in {1,N-1}, Q=0, 101/125<gamma<27/25,

  with own native5, own native>=7 or dotted.

The dotted upper bound uses m>=5 and
sigma^2/delta>=4*(27/25)-2>(761/500)^2 at gamma=27/25,
exactly as in the pinned general-N proof. The near4 role with
helper-side m>=8 is strong. At m7, q_7 is in (3,25/8) for BOTH
N20 and N21, so the sharper ratio F(q)/q<3/2 makes that role
strong too. A weak near4 thus has its actual mate at4 or5,
or at the reflected indices. The mate's normalized height is
greater than sqrt3/2>43/50, by the N>=15 actual-mate lemma.

## 3. Adjacent weak literal far-helper caps are still impossible

This step is proved before removing near4 roles. Consider two adjacent
literal far-helper caps at i,j=i+1. Put f=101/125, g=27/25,
A=gamma_i, B=gamma_j, and write their positive time coefficients
as a_i,a_j. Keeping the helper terms gives the useful smaller-rank
time bound

    a_i^2 > (4/5)^2*(15/2)+(4/5)*(1/2)+19/40-1
            =187/40>4,

and the same for j. This replaces the older bound that discarded
the helper terms and required D_i>=85/11.

The weighted AM--GM expansions from the pinned general-N proof are
identities at N20/21 as well. For the actual mutual product R they
give

    same helper:
      R>=-AB(1-1/N) +/- (A-B)/N +(N-1)/(2N)
                           +(N+1)(A/B+B/A)/(4N)>1-g^2>-1/2;
    outward opposite helpers:
      R>=1/2-AB+(A-1)(B-1)/N;
    reversed opposite helpers:
      R>=1/2-AB+(A+1)(B+1)/N.

The two opposite-helper lower bounds exceed
1/2-g^2-(1-f)(g-1)/20>-7/10>-1/sqrt2.
The native gap consequently gives R=0 for a common helper and
R in {0,-1/2} for opposite helpers.

For outward opposite helpers R=0 is impossible, because eta_i,eta_j<4/5
and the exact determinant comparison is

    (D_i-4/5)(D_j-4/5)-C_ij^2
       =(2/5)C_ij-24/25+8/(5N)>0,

already using C_ij>=7. For the reversed assignment put l=i,
s=N-j, so l+s=N-1 and l,s>=5. The exact determinant polynomial is

    NF=4ls A^2 B^2-(l+1)(s+1)(A^2+B^2)-2AB
                        -(l+1)A-(s+1)B+(N+2)/4.

Its derivative in either A or B is positive on [4/5,g]: for example
the positive bracket is at least (39/25)ls-N, and ls>=70,
N<=21 makes its derivative positive even after subtracting
2g+(l+1). If one gamma is at least9/10, concavity in l gives

    NF >= (509N-9974)/500 >=206/500>0.

The value follows by substituting l5 or s5 and (A,B)=(4/5,9/10)
or its reverse. It remains valid at N20, not just N>=22.
If both gammas equal phi/2, the pinned literal far5 distance-one
prohibition applies because it is proved for every N>=15.
Hence opposite helpers have R=-1/2 in every possible case.

The actual two-cap replacement and its full displacement
Z=W_i omega_i+W_j omega_j, W_i,W_j>0, now exist exactly as in
the pinned mixed-pressure source. Append all old spectator normals;
they form an orthogonal positive block and the pressure points have
zero spectator displacement. They remain feasible against every
unknown actual wall. Put W=W_i+W_j and
E0=W_i sqrt(D_i)+W_j sqrt(D_j). The unchanged two-by-two Schur
calculation gives

    sqrt(1+Dcal)/W < 7/16     for a common helper,
    sqrt(1+Dcal)/W < 3/10     for opposite helpers.

For completeness, the first bound uses only l>=5, a_i+a_j>4,
and W>(4/5)l(a_i+a_j). The second uses delta_l,delta_s in
(1/2,3/5], 5/6<lambda_i/lambda_j<6/5,
Dcal/(lambda_i+lambda_j)^2<=3/55 and lambda_i+lambda_j>40.
These are the same explicit inverse calculations, with both
opposite-helper orientations retained.

The other displacement estimate needs a slightly weaker constant
at the smaller rank. Rationalization gives

    E0-||Z|| < W/[2(sqrt(15/2)+sqrt7)]<W/10.

The last inequality follows from sqrt(15/2)>27/10 and sqrt7>13/5,
whose sum exceeds5. Thus the actual mixed-point inequality bounds
the weighted average of the two losses by

    7/16+1/10=43/80<61/100     for a common helper,
    3/10+1/10=2/5<61/100      for opposite helpers.

Two weak literal caps would give both losses at least61/100, a
contradiction. This proves their nonadjacency with all spectators.

## 4. Removing near4 under the contrary energy bound

Suppose a weak near4 cap occurs on the left. Its actual mate is at
j<=5 with normalized demand greater than43/50. If a weak near4 cap
occurs on the right, its reflected mate k>=N-5 has demand greater
than43/50. Otherwise both indices N-6,N-5 are either strong or
literal far-helper caps, so Section3 supplies a strong one k.
Its normalized demand is greater than

    1-(61/100)/(27/10)>77/100.

In either case k>=N-6, and the dual-direction correlation obeys

    rho=sqrt(j(N-k)/((N-j)k))
         <=sqrt(5*6/(15*14))=1/sqrt7<1/2.

Consequently their two normalized projections sum to more than
43/50+77/100=163/100, while the squared norm of the sum direction
is less than3. Cauchy and e<=6/7 are impossible, since

    (163/100)^2-3*(6/7)=5983/70000>0.

Reversing the normal numbering treats a right weak near4 cap.
All remaining weak central directions are literal far-helper caps;
none are adjacent.

## 5. Two strong directions spend more than six sevenths

There is a strong i in {5,6} and a strong j in {N-6,N-5}.
Their first and last coordinate intervals are disjoint. The exact
height energy and Cauchy on those intervals give

    e >= h_i^2/(2i)+h_j^2/(2(N-j)).

At N20 the endpoint floors are already

    q5-61/100 > 21/10,
    q6-61/100 > 57/25,

and increasing N improves them. Positive squaring verifies the
bounds: 15/2>(271/100)^2 and 42/5>(289/100)^2. By reflection the
same estimates apply to j. Each displayed tail-energy term exceeds

    min((21/10)^2/10, (57/25)^2/12)=1083/2500.

It follows that

    e>1083/1250>6/7,

with final difference81/8750>0. This contradicts the contrary
assumption and proves the theorem.

## 6. Restricted application and verification boundary

For a hypothetical compact dimension23 farthest forest, the theorem
excludes A19 plus rank4 spectators and A20 plus rank3 spectators.
Their ordinary ranks cost rank/14, with at most the single permitted
short-A1 discard or instead an isolated rank-two7 block of cost>9/100.
The weakest short-A1 alternative for A20 already gives

    6/7+2/14=1,

strictly exceeded. The rank-two alternative gives
6/7+1/14+9/100>1. The A19 cases are stronger. This use requires
the separately reviewed universal exception/ordinary-cost theorem.
It does not exclude the remaining A16 or B/D frontier by itself.

The new ingredients are the smaller-rank weak-inventory seals, the
N20/21 adjacent determinant bound, removal of near4 at the stated
contrary projection bound, and the disjoint endpoint energy. All
actual-point, classification and spectator tools retain their
documented contributions. This note awaits source-bound full root/
outside reviews. No checker output or new Lean coverage is claimed.
