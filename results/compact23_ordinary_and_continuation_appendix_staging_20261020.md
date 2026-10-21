# Compact23 companion appendix: ordinary costs, small debits and continuation methods

20 October2026. Staging prose for the existing dimension23 manuscript,
prepared by resource_cleanup. No primary TeX, accepted source, manifest,
catalogue or frozen packet is changed; no computation, import, Lean or
TeX job was run. Pending dimension22 candidates are NOT inputs here.

The proposed ordinary-cost appendix below supplies the low-rank proofs
behind Proposition ordinary, with notation shared by the manuscript.
The final two sections group Q and the large caps by proof method and
identify the exact finite calculations still requiring expansion. They
do NOT pretend that a short source map makes every companion derivation
part of a standalone paper. This is one supplement to the existing paper,
not a seventh paper or a new mathematical result.

## A.1. One common mechanism; arbitrary spectators retained

Use the manuscript's actual farthest vector z=E(p-x), ||x||^2<1,
old unit Gram B, inverse C, d_i=C_ii, omega_i=sum C_ij*n_j and
u_i=omega_i/sqrt(d_i). For an actual outside wall c=A p-v put
beta_i=-c.n_i>=0. Its component cost is q=beta_T^T C_T beta_T,
its total foreign cost is Q, and A^2=q+Q-1>0. Write
v_i=(C beta)_i/sqrt(d_i) and r_i=tanh L_i.

We use the manuscript's proved actual-wall alternatives, not a formal
completion: v_i/A>1 makes c the actual opposite of i; otherwise a
distinct opposite forces r_i>=v_i/A. Equality gives an ideal endpoint
and is impossible compactly. Distinct actual opposites therefore obey
r_i*r_j>=K_ij=u_i.u_j. If a component has energy at most e, then

    r_i*r_j<=4s/(1+s)^2, s=e(1+K_ij)/2.                 (A1)

Indeed the half-tangent sum is at most sqrt(2e(1+K_ij)), and
log(2t/(1+t^2)) is concave on(0,1). A strict reverse comparison in
(A1) forces the SAME actual opposite. This conversion must retain
the distinction between full and half tangents.

If c closes an edge in T, another whole component is retained in that
actual finite replacement. It has at most one c contact: otherwise
the retained diagram has a nonpositive cycle. Each nonzero foreign
component cost is at least1/4. In a nontrivial component its inverse
diagonal is at least4/3 by restricting to one adjacent3 bond; thus a
cost strictly between1/4 and1/3 does not occur. Most importantly,

    every contacted foreign component has cost strictly below A^2. (A2)

Otherwise its contacted derivative has square at least its cost,
so the reciprocal shadow makes c close an orthogonal whole-component
direction too, or gives a null endpoint. Both are forbidden by the
manuscript's orthogonal-opposite lemma. Consequently A^2<1/4 implies
Q=0. This universal isolation argument replaces spectator catalogues.

When one wall closes ALL directions of a component, put h=A^2 and
R_e^2=4e/(1+e)^2. Then

    q-1<=h<=R_e^2*min_i v_i^2.                          (A3)

Contacts whose connecting old path omits another old node are impossible:
delete that omitted direction, also closed by c; its retained actual
finite diagram contains the cycle. In a path component multiple
contacts may therefore occur only at its two endpoints. In D4 there
can be only one contact. A single-contact own cost q=1 gives a positive
affine kernel in the component+c principal set, forbidden compactly
even when further spectators are present.

All inverse-Gram tables below follow by inversion of the unit path/star
Gram with the stated bond labels. Thus they are exact displayed
algebraic identities, not an unmentioned numerical census. For a single
contact define sigma=min_i(v_i^2/q); (A3) gives
q<=1/(1-R_e^2*sigma).

## A.2. Short A1 directions and the complete debit rule

Call L<asinh(1/sqrt2) short. If its closing wall has released square
b and squared height h, then

    h=b+Q-1, h<b/3, Q<1-2b/3.

There cannot be two contacted foreign components: Q>=1/2 would imply
h<(1-Q)/2<=1/4, violating(A2). With one, Q<h=b+Q-1 forces b>1,
and Q<1/3. Hence the contact is precisely A1_3 of cost1/4. Without
foreign contact, h=b-1>0 likewise makes the release dotted. The ONLY
short A1 types are therefore

    I: isolated dotted, 0<h<1/2;
    J: dotted plus one A1_3, initially1/4<h<3/8.

For J, strict interiority improves the lower bound to h>1/3. The
heavy edge contributes more than A/2 to x.v. On the light A1 the
shadow contributes at least A*(1-sqrt(1-1/(4h))); if h<=1/3 this is
at least A/2. Their sum contradicts x.v<A. It follows that the heavy
direction has r^2>4/13 and half-tangent square>9/100, while the light
direction has r^2>2/3>16/25 and half-tangent square>1/4. Thus J funds
its two A1 factors above9/100+1/4>2/10.

There is at most ONE short A1 direction. Two such distinct walls have
height product<1/2 and nonnegative spatial product. The actual native
gap therefore makes them orthogonal. Two I/J supports cannot realize
this: disjoint supports give a negative product; a heavy overlap gives
spatial product>1/2 (indeed>1 for heavy/heavy). For two J walls the
sole possible light/light overlap gives1/4, LESS than their height
product because both squared heights exceed1/3. Every case is nonzero.

A nonshort A1 has energy at least

    tanh^2(asinh(1/sqrt2)/2)=5-2sqrt6>1/10,

using sqrt6<49/20. Thus only I can cause an unfunded A1 debit; it is
always safe to discard just that one rank. Alternatively a J pair
may be counted at the stronger ordinary two-rank budget.

The only possible nonordinary rank-two blocks will be proved in the
next section to be OLD A2/NEW7 and OLD I2(7)/NEW3. Their squared
heights are respectively in(0,1/12) and(0,1/3), and each pays>.09.
Their spatial supports are wholly in their own old blocks. Two such
defects in different blocks, or one defect and any I/J, have actual
wall product -sqrt(h1*h2) in(-1/2,0), impossible. Therefore at most
one short-A1 debit OR one isolated rank-two defect can occur, never both.
This argument is independent of all other finite component types.

## A.3. A2 and the entire dihedral range

For A2 put e2=319/2000 and sinh^2 L0=31/50. If both actual lengths
are at least L0, their half-tangents exceed173/500; correlation1/2
gives energy>4*173^2/(3*500^2)>e2. A shorter direction cannot have
two positive own contacts:3/3 is affine, while3/4 or larger gives
sinh^2 L>=4sqrt2/9>31/50, monotonically increasing with either
contact and Q. If only the other node is contacted, classify that
shorter contacted direction first. For its single contact square b,

    h=4b/3+Q-1, 0<h<(124/243)b,
    Q<1-(200/243)b.                                    (A4)

Every positive foreign cost must be<h. Label3 has h<1/4 and cannot
give positive height after isolation. Label4 has h<62/243; every
foreign cost must then equal1/4. Only Q=2/4 could make h positive
within the bound, but it gives h=1/6<1/4, impossible. Label5 has Q<1/2,
so at most one foreign component, whose cost exceeds h since the own
cost phi^2/3<1. Label6 has an affine3,6 principal.

For native>=7, b>.81 gives Q<1/3: Q=0 or one A1_3. The latter can
be short only for b<729/800, exactly m7,8,9,10; for m>=11,
cos^2(pi/11)>45/49>729/800 by pi<22/7. For dotted b>1, (A4) gives
Q<1/4, so Q0. The surviving possibilities are therefore isolated
finite>=7/dotted, or A2_m+A1_3 for m7..10.

For isolated finite m, the same wall closes the mate since its
reciprocal ratio square is b/(4b-3)>1. At m>=8 the ACTUAL mate has
r^2>=6sqrt2-8>400/841, hence half-tangent>2/5 and block energy>4/25>e2.
For isolated dotted b, (A4) gives b<243/200; the DISTINCT mate shadow
has r^2>=b/(4b-3)>81/124>16/25 and energy>1/4. For mixed m7..10,
the mate cannot be closed by this wall (its white cost b+1/4>1);
the shadow gives r^2>=b/(4b-9/4)>11/17>16/25, again energy>1/4.
Only isolated NEW7 at OLD A2 remains. Its h=4b/3-1 lies in(0,1/12).
The cubic root t=cos(2pi/7) satisfies623/1000<t<5/8, verified in
8t^3+4t^2-4t-1, increasing on[1/2,1]. Thus b>1623/2000 and the
actual mate has r^2=4-3/b>164/541>3600/11881; its energy exceeds9/100.

Now let the OLD block be I2(m), m>=4. Set k=cos(pi/m), d=(1-k^2)^-1,
so C=d[[1,k],[k,1]]. Under energy<=1/7, (A1) forces one opposite
for both directions. If both contacts are positive, the actual edge
has cosh L>=cot(pi/(2m))/sqrt3>11/8; its individual energy exceeds
3/19>1/7. Thus there is one contact beta, own cost q=d*beta^2, and

    q*(1-7k^2/16)<=1.                                 (A5)

The minimum beta^2=1/4 makes the left side(9d+7)/64. For m>=8,
d>=4+2sqrt2>19/3, so no such contact exists. This proves a finite
old-label reduction, rather than assuming one.

| Old m | Exact consequence of(A5) | Conclusion |
| --- | --- | --- |
| 4 | beta^2<=16/25<cos^2(pi/5); labels3/4; h<=7/32<1/4 | isolation; q<=1 |
| 5 | d>25/9 excludes beta>=1/sqrt2; only3; h<7/32 | isolation; q<1 |
| 6 | only3, q1 | affine kernel |
| 7 | d>5 excludes beta>=1/sqrt2; Q<3/16 | Q0; isolated NEW3 |

The last case has h=d/4-1 in(0,1/3). Its actual mate has
cosh L=cot(pi/7)/sqrt3>109/91, using the same cubic bracket for t;
the final squared comparison has positive numerator2652. Its energy
is therefore>9/100. This is OLD I2(7)/NEW3, distinct from OLD A2/NEW7.
The native and dotted beta ranges were all included in(A5).

## A.4. Rank-three and rank-four costs in one calculation

At e1/4, (A1) forces a common opposite along adjacent correlations
>=1/sqrt3; the cleared lowest comparison is74/sqrt3>39. At e1/3
the analogous sufficient correlation is1/sqrt2, with51/sqrt2>34.
For A4 at e2/7 the adjacent correlations are sqrt(3/8),2/3,sqrt(3/8);
the cleared comparison is k^3+16k^2+36k-28>0, whose lowest value is
positive because291sqrt6>704. For H4 use e11/36<1/3. Direct inversion
of each indicated finite unit Gram shows its adjacent correlations
meet these thresholds. Hence one opposite closes the ENTIRE component.

Let the path bonds be A3:(3,3), B3:(4,3), H3:(5,3), A4:(3,3,3),
B4:(4,3,3), F4:(3,4,3), H4:(5,3,3); D4 is the simple three-arm star.
The following endpoint-contact tests are exact. Two endpoint3 contacts
in A3/A4 are affine cycles; otherwise their next minimum has one4.
For the other rows both contacts may already be3. Each expression
is q-R_e^2*v_j^2, which(A3) says is at most1.

| Type | Test direction j | Lower bound for two endpoint contacts |
| --- | --- | --- |
| A3 | middle | 177/200+9sqrt2/100>1 |
| A4 | second | (230+15sqrt2)/243>1 |
| B3 | middle | (77+18sqrt2)/100>1 |
| B4 | third | 15/16+sqrt2/8>1 |
| F4 | second | 15/16+sqrt2/4>1 |
| H3 | middle | (72phi+95)/100>1 |

D4 has no possible pair of contacts by the omitted-node cycle argument.
For H4 two contacts could only be its endpoints, but deleting old4
then retains a path at least3,5,3 of determinant(5-4phi)/16<0.
Thus all these components have a SINGLE contact.

For the crystallographic and A cases the complete single-contact
calculation is the following table. Repeated entries are retained
in path order; D4 uses leaf,center. R_e^2 is16/25 at e1/4,56/81 at
A4 e2/7, and3/4 at e1/3. The upper q follows directly from(A3).

| Type | d_j | sigma_j=min(v_i^2/q) | Upper q |
| --- | --- | --- | --- |
| A3 | 3/2,2,3/2 | 1/9,1/3,1/9 | 225/209,75/59,225/209 |
| A4 | 8/5,12/5,12/5,8/5 | 1/16,1/6,1/6,1/16 | 162/155,243/215,243/215,162/155 |
| B3 | 3,4,2 | 1/3,1/2,1/3 | 75/59,25/17,75/59 |
| D4 | 2,4 | 1/4,1/2 | 16/13,8/5 |
| B4 | 4,6,4,2 | 1/4,1/3,1/2,1/4 | 16/13,4/3,8/5,16/13 |
| F4 | 4,12,12,4 | 1/2,2/3,2/3,1/2 | 8/5,2,2,8/5 |

Dividing each q upper by d_j bounds the native contact square.
The resulting possibilities all have q<=1 (or none exists). To
make this finite check explicit: A3 endpoints allow3/4/5, middle3/4;
A4 endpoints3/4, inner3; B3 nodes allow3;3;3/4; D4 leaf3/4,center3;
B4 allows3;none;3;3/4; F4 endpoints3,innernone. Every q1 is affine.
For every other case(A3) bounds h<1/4: the largest relevant bounds
are16/75,7/155,28/405,3/16. Isolation then gives h=q-1<0. This proves
A3,B3>1/4, A4>2/7 and B4,D4,F4>1/3.

For H3, at minimum contact beta1/2 the q-R^2*v_i^2 tests for contacted
nodes1,2, with test direction3, are respectively
(199phi+263)/500 and(61phi+77)/125, both>1. For contacted node3,
beta>=1/sqrt2 gives test at1 at least(27phi+118)/150>1. Only3@3
remains, with q=(phi+2)/4<1 and h<=4phi^4/75<phi/4. No spectator
gives negative height; one gives Q>h; two give Q>=1/2 and therefore
h>=q-1/2=phi/4. All are impossible, proving H3>1/4.

For H4 set R^2=1584/2209. At beta1/2 the four single-contact tests
(contact j, test i) have values

    (1,4): (4938phi+3688)/2209;
    (2,4): (7209phi+5334)/2209;
    (3,1): (2500phi+2667)/2209;
    (4,1): (625phi+1219)/2209.

All exceed1, the last from phi>198/125. Larger beta scales these
positive tests upward. Thus H4>11/36 universally. All proofs in this
section retained arbitrary spectators using(A2), not an assumed seal.

## A.5. A5 and A6: complete finite-white contact analysis

Use the path inverse C_ij=2min(i,j)(N-max(i,j))/N, N6 or7.
For A5, distinct2/3 and3/4 pairs each demand sum>109/100. The vector
u2+2u3+u4 has norm^2=7+4sqrt2<443/35, and
(109/50)^2-(3/8)*(443/35)=209/35000>0. Thus under energy<=3/8
there is, up to reversal, a common2/3 opposite.

For A6 the pairs2/3,3/4,4/5 have distinct-wall sums>2249/2000,
11547/10000,2249/2000. The vector u2+2u3+2u4+u5 has norm^2
<33781/1250. The demand is34037/10000, with Cauchy gap
2181583/700000000 at e3/7. Hence a common2/3 or3/4 occurs.
These rational pair-sum seals follow by substitution in
F(s/2)^2<K; no approximate optimization enters.

A common consecutive pair permits one contact on either retained
path, so at most two contacts. The COMPLETE double-position lists are:

| Block/pair | Positions | Necessary finite-white labels |
| --- | --- | --- |
| A5,2/3 | 1,3 | 3/3 or4@1+3@3 |
| A5,2/3 | 1,4;1,5;2,4;2,5 | both3 |
| A5,2/3 | 2,3 | both in3/4/5 |
| A6,2/3 | 1,3 | 3/3 or4@1+3@3 |
| A6,2/3 | 1,4;1,5;1,6;2,4;2,5;2,6 | both3 |
| A6,2/3 | 2,3 | 3/4/5 at2,3/4 at3 |
| A6,3/4 | 3,4 | both in3/4/5 |
| A6,3/4 | other8 left/right pairs | both3 |

To verify exhaustiveness, delete one closed direction: every listed
nonadjacent double survives in a connected rank5/6 finite path or
D/E6 tree. A4 can occur only at a path end; D/E are all3. The sole
4@1 row is precisely the permitted terminal4. All3 doubles give
an affine cycle. Adjacent non3/3 doubles have an actual endpoint
cosh at least(3+sqrt2)/2>11/5 in A5. In A6 their minimal formulas are
(5+4sqrt2)/sqrt7, (4sqrt2+6)/sqrt14, or(12+9sqrt2)/sqrt42,
each>5/2. These exceed the individual-energy thresholds for3/8,3/7.
The A5 4@1+3@3 row has Q0 and cosh L3=(3+sqrt2)/2. The A6 row
has Q0 (white cost9/10) and is handled by its actual foot below.

Every single contact is3/4/5; a6 with an old3 is nonpositive in one
retained endpoint. For single3 write h=q+Q-1. The complete q,Q bounds
are:

| Block/pair | Contact j | q | Upper Q |
| --- | --- | --- | --- |
| A5,2/3 | 1,2,3,4,5 | 5/12,2/3,3/4,2/3,5/12 | 2/3,2/3,5/8,1/2,5/8 |
| A6,2/3 | 1,2,3,4,5,6 | 3/7,5/7,6/7,6/7,5/7,3/7 | 2/3,2/3,3/5,2/5,2/5,3/5 |
| A6,3/4 | 1/6;2/5;3/4 | 3/7;5/7;6/7 | 5/8;1/2;5/8 |

Every end or small-Q row has h<1/4, or allows just one foreign
component whose Q>h. In the other A5 rows two components would each
exceed1/3 (j2) or1/4 (j3), forcing Q>=2/3. In A6 j2 the same cost
gap applies. For the remaining A6 j3/4 rows two foreign components
must have first bonds3, as Q<5/8; c, its internal contacted old node,
the two old neighbors and those two foreign nodes form exact affine
two-fork D5. Thus all single3 cases disappear. A single4 internally
is the affine3,3,4 star. At an end the white Q bounds are <=1/3
and h<1/4, contradicting(A2) or isolation.

A single5 can attach only at a retained path end; its H3/H4 white
cost is>3/4, so Q0. The following actual directions cover ALL cases:

| Block/pair | Contact | Actual direction | cosh L |
| --- | --- | --- | --- |
| A5,2/3 | 1 or5 | 4 or2 | phi^3/sqrt3 |
| A5,2/3 | 2 | 3 | phi^2 |
| A5,2/3 | 3 | 2 | sqrt3*phi^3 |
| A6,2/3 | 1 | 4 | 3phi^3/sqrt21 |
| A6,2/3 | 2 | 3 | 4phi^2/sqrt14 |
| A6,3/4 | 1 or6 | 4 or3 | 3phi^3/sqrt21 |
| A6,3/4 | 3 or4 | 4 or3 | 9phi^3/sqrt21 |

The extra endpoints use positive full replacement of the ALREADY
actual wall; no new facet is inferred from a formal Gram. The A5
entries exceed11/5 and the A6 entries5/2, by phi>8/5,
sqrt21<23/5 and sqrt14<15/4. Finally the A6 double4@1+3@3 gives

    xi=c+(sqrt2-1/2)n1+(sqrt2-1)n2+(sqrt2/2)n3,
    xi.c=xi.n1=xi.n3=0, xi.n2=(sqrt2-3)/4<0,
    xi^2=5/4-sqrt2=-k, A^2=(5+4sqrt2)/7.

It is a positive actual future point against every wall. Its pure-A
pressure has energy at least(A-sqrtk)/(A+sqrtk)>1/2 because
A^2-9k=(335-236sqrt2)/28>0. This completes A5>3/8 and A6>3/7.

## A.6. B5 through B8: one forced pair and every short-tail bridge

Use u_i=prefix_i/sqrt i, omega_i=sqrt2*prefix_i for i<r,
omega_r=prefix_r. At energy3/7 the pair2/3 is common for r5/6;
at1/2,4/7 the pair6/7 is common for r7/8. Substitution in(A1)
gives cleared signs127sqrt(2/3)-100, (69sqrt(6/7)-52)/7 and
(199sqrt(6/7)-176)/7, all positive.

For B5/B6 common2/3, deleting2 retains A1+B3/B4. Its right finite
contact is3 at long end3, with the sole B5 exception3 at short5,
making F4. Deleting3 retains left A2, so at most one left3/4/5.
For long3@3, zero own2 is affine and own3 is an affine triangle;
own>=4 gives cosh L3>=2+3/sqrt2>3. With no right contact, the
complete left choices are3/4/5 at1 or2:3@2 and both4 choices are
affine;3@1 has h=Q-1/2<1/6 and forces a forbidden foreign closure.
Both5 choices have Q0;5@1 actually closes4 with cosh L4=phi^3>3,
and5@2 gives cosh^2 L3=2phi^2/(2-phi)>9. The B5 short5 contact3
seals left1 and all spectators by finite F4; own2=0 is affine, and
the other actual endpoint forces own2=3, giving cosh L2=2+sqrt2.
Thus all roles exceed the necessary individual threshold11/5.

For common6/7 in B7/B8, no right contact leaves a single left3 or
endpoint4. All but far3@1 are affine fork/double4/star principals.
A far3@1 with arbitrary spectators is also impossible: a minimal
nonpositive principal must contain c and7, since deletion7 is actual
finite; their only old connecting path c,1,...,7 has8 vertices.
Spectators attach only through c and do not shorten it. This violates
the Lanner order bound5. Thus the no-right branch is FULLY excluded.

In B7 the retained right node is an A1 with finite contact beta.
The complete possibilities give:

| beta/own6 | Other left contacts | Actual pressure or rejection |
| --- | --- | --- |
| beta>=cos(pi/5), any own | none by finite rank7 classification | cosh L6>=sqrt6*cot(pi/5)>3 |
| beta4, own0 | any allowed | affine4,4 |
| beta4, own>0 | none | cosh L6>=sqrt24 |
| beta3, own>0 | none | cosh L6>=2+sqrt2 |
| beta3, own0 | none or far1 | exact affine F4 |
| beta3, own0 | one j2..5 | cosh L6>=1+3/sqrt2>3 |

For the last row the FULL white data are
Delta<=3/4-j(6-j)/12<=1/3 and
v6>=sqrt12*(j/12+1/(2sqrt2)). They retain every foreign Q and
all E7 extension positions. Larger contacts only strengthen the bound.

In B8 the retained right B2 has exactly one simple3 at either end.
At long end7, zero own6 is affine; own3 is an affine triangle; own4
gives cosh L7>=6+7/sqrt2>7. At short end8, own6>0 gives
cosh L6>=sqrt6+sqrt3>11/3. For own6=0, the remaining left-white
budget forces a3 at1 or5 if any contact exists: none/1 are affine F4,
whereas3@5 is an ALLOWED bridge, with

    Delta<=1/12, v6=sqrt12*(5/12+1/(2sqrt2)),
    cosh L6>=5+3sqrt2>9.

No false F4 sealing is used for this bridge. These lists exhaust both
actual positive replacements, including every zero-own and foreign role.
Their edge thresholds prove B5/B6>3/7, B7>1/2 and B8>4/7.

## A.7. Extended D4, the remaining D7 case and the D8 analytic closure

An old central D4 star S has center o and three simple neighbors.
Suppose one arm extends and a single actual wall closes ALL four
directions of S. A minimal nonpositive principal in the old vertex+c
must contain c and every node of S: deleting any S-node is an actual
finite replacement. Compactness removes affine cases and the Lanner
order bound forces this principal to be precisely S+c.

Its complete Lanner classification is elementary. Two c contacts give
a nonpositive proper triangle or cycle. A sole center3 gives affine
D4; center4 or larger has a nonpositive proper star. A sole leaf contact
can only be3/4/5 from its finite proper path;3 is finite D5,4 affine,
so only leaf5 remains. If that leaf extends, c and the extension give
the nonpositive5,3,3 star; if another arm extends, they give the
nonpositive5,3,3,3 path. Each omits an S-node and survives its supposed
finite replacement, impossible. Extra negative contacts cannot repair
the positive nonpositive quadratic witness. Thus the four directions
cannot share one opposite, with arbitrary spectators.

For arms(1,2,L), L2,3,4 in E6,E7,E8, the raw norm2 central Schur
residual is delta=2-[1/2+2/3+L/(L+1)]. A neighbor on arm length m
has squared central correlation t/(delta+t), t=m/(m+1). The minima
are sqrt3/2,sqrt(6/7),sqrt15/4. At energies9/20,5/9,2/3, even the
coarse pair upper4e/(1+e)^2 is below each. Thus all four opposites
would coincide, a contradiction. For D5,D6 the minimum central
correlations are sqrt(3/5),sqrt(2/3); (A1) at e5/14,3/7 has cleared
values544sqrt(3/5)-362 and127sqrt(2/3)-100, both positive. This
proves D5>5/14,D6>3/7,E6>9/20,E7>5/9,E8>2/3.

For D7 let o=r-2,a=r-3 and fork leaves b,e. Under energy<=1/2,
(A1) at correlation sqrt(4/5) forces one opposite for o,a. The full
finite-white classification has exactly three complete profiles:

    4@a+3@o; isolated5@b; 3@a+4@b, up to swapping b,e.

Here is the exhaustive argument. Deletion a retains A3=b-o-e, so
the right contact is none, middle3, or endpoint3/4/5. Deletion o
retains one long A path and permits at most one left contact.
Middle3 forces4@a (own0/3 affine); endpoint5 seals H4 against all
white extensions, and a contact to the removed a instead gives the
retained nonpositive3,3,3,5 path. Thus it forbids every left/spectator
extension. Endpoint4 forces3@a (otherwise affine
terminal4 fork). Endpoint3 forbids every left contact by cycles or
retained F4 witnesses; a foreign3 gives affine E7, foreign4 gives an
affine terminal4 fork, larger foreign bonds are already nonfinite.
It is then isolated3 at a fork: for r7 its height is negative, and
for r>=8 it contains affine E8. With no right contact, internal-left3
and endpoint4 are affine; far-left3 requires a connected Lanner with
at least6 vertices along c to o, impossible. This includes ALL foreign
contacts; the finite H4/B4 end cannot be extended at c.

The first/third profiles have cosh L_o=a+(r-2)/sqrt2>7. The isolated5
profile has r_a^2=1-(9-4sqrt5)/(r-3); at r7 its half-tangent square
is4/sqrt5-1>3/4. Thus every common o,a profile costs>3/4, contradicting
the assumed1/2 and proving D7>1/2.

For D8, common4/5 has exactly4@4+3@5, and common5/6 is covered by
the preceding center/neighbor list. To prove the first assertion,
deleting4 permits one simple3 at a right D4 leaf, and deleting5
permits one left3 or endpoint4. Long-leaf3 forces4@4 and seals
foreign contacts by the retained3,4,3,3 path. A fork3 and any left3
give an affine cycle; no left gives affine E8; left4 is retained
nonpositive. With no right contact, only far3@1 survives the affine
tests, and the Lanner path to5 already has6 vertices. Thus both common
pairs force an actual half-tangent square>3/4. At energy<=4/7 the
actual4/5 and5/6 opposites are therefore distinct.

For completeness the remaining three-direction analytic contradiction
is recorded explicitly. Put q=sqrt(4/5), k=sqrt(5/6), lambda=(1+k)/2,
alpha=k*sqrt(1-k)/2, e=4/7 and P(t)=t/(1+sqrt(1-t))^2. Write
x=a*w+b*v+z*t0+y where w,v are the normalized sum/difference of
u5,u6 and u4=q*u5+t0/sqrt5. The two full-tangent inequalities give

    h4*h5>=P(q),
    T=(1+alpha)b^2+lambda*z^2<=delta=lambda*e-P(k)<1/400.

The latter follows by the tangent line to4p/k-(1+p)^2 at its first
zero p=P(k), using p=h5*h6=lambda*a^2-(1-lambda)b^2.
Retaining the negative quadratic terms in h4*h5 gives

    h4*h5<=q*lambda*e+A|b|+B|z|-q*k*b^2-q*lambda*z^2,
    A^2=8/105, B^2=4/35.

Weighted Cauchy bounds the linear part by sqrt(C*T), where
C=(8/105)/(1+alpha)+(4/35)/lambda<47/250. Also q*k/(1+alpha)>.71,
q>.71 and q*lambda*e<.489. Hence

    h4*h5<.489+sqrt((47/250)T)-.71T
            <101781/200000<.509<P(q),

a contradiction. All constants are exact rational seals: .8944<q<.8945,
.9128<k<.9129, .29<sqrt(1-k)<.2953, .13<alpha<.14,
.95<lambda<.95645. The tight delta inequality reduces to
91280000*70000-167780209*38083=26300653>0; the square-root bound in
the last line is sqrt(47/100000)<271/12500. Thus D8>4/7 is proved,
without turning a higher-rank energy rate into a rank-proportional one.

Together A.1--A.7 supply the full ordinary table and its precise small
debit exceptions. The last two sections below concern different inputs.

## B. The Q input: finite continuation first, actual feet second

This section organizes the B/D>9/14 proof as two mathematical lemmas
followed by one common numerical implication. Its complete low-index
pressure algebra remains in Q1/Q2 and their named sublemmas; the end of
this note identifies which parts have not been fully transcribed here.
The six rational dual coefficients already displayed in the manuscript
do not by themselves prove distinctness of the three actual pairs.

### B.1. Complete common-pair continuation lists

Use the standard B/D unit normals and u_i=prefix_i/sqrt i for i<=7.
Under the assumed actual component energy<=9/14, common4/5 or6/7
internal pairs have exactly4@i+3@(i+1), ALL other old contacts0, and
cosh L_(i+1)=i+(i+1)/sqrt2>7. Here the right retained B/D rank is>=4.
The proof is the internal-pair version of A.7: finite extension allows
one right long3; zero-own is affine; own3 is a triangle; own4 seals
all spectator contacts by a retained3,4,3,3 path. A fork contact in D
gives an affine E8 or a cycle, and no right contact gives an overlong
Lanner path. These claims use actual finite replacements, not the
false assertion that every nonpositive principal is forbidden.

At r9 there are two genuine boundary cases for6/7. In D9 use the
center/neighbor profiles of A.7, each forcing energy>3/4. In B9 the
right B3 permits the short-end3 extension to F4. It seals every left
and foreign WHITE contact, but the removed own6 must still be checked.
Own0 is affine F4; own4 is the affine3,3,4 star; the remaining own3
has cosh L6=(1+1/sqrt2)*sqrt12>11/2, giving energy>9/13. This proves
distinct6/7 in every B/D rank>=9 without dropping the r9 F4 role.

For common1/2, deletion1 retains B/D rank>=8: its only nonzero finite
attachment is3 at long leaf2. If present, own0 is affine, own3 an
affine triangle, and every own>=4 has

    cosh L2=(1+gamma)/sqrt(1-gamma^2-Q)>=1+sqrt2.       (B1)

With no right contact, the complete profile is a single finite
gamma at1 and arbitrary foreign cost subject to

    A^2=2gamma^2+Q-1>0, gamma^2+Q<1.

Own0 has zero derivative; own3 has height<1/4 and forces a foreign
closure; own4 and6 are affine. For own5 and Q>0, Q>=1/4 gives
cosh L2>=phi^2>1+sqrt2. With Q0 it is isolated5. For own>=7,
the endpoint condition forces Q0; own>=8 satisfies the threshold(B1).
The ONLY last role is isolated native7 at1. This completes every
1/2 foreign branch before using the pressure lemma below.

For common2/3, deletion2 retains B_(r-2) or D_(r-2), and deletion3
retains left A2. The latter allows precisely one3/4/5 at1 or2.
The right long3@3 forces a positive contact at2: own0 is affine;
own3 affine triangle; own4 gives4@2+3@3; own5 has Q0 from white cost
phi^2/3>3/4 and cosh L3=(2+5phi)/sqrt2>7. With no right contact,
the exhaustive left list is

    3@1: A^2=Q-1/2<1/6, hence forbidden foreign closure;
    4@1: terminal4 affine;
    5@1: isolated5 after Q0 from its actual A2 white cost;
    3@2: affine two-fork D or terminal4 B fork;
    4@2: affine3,3,4 star;
    5@2: isolated5 after the same Q0.

For D only, a right fork3 can occur at the r9 boundary: the retained
D7+c is finite E8. This DOES forbid extra contacts to old1 and any
spectator, since they would enlarge that retained connected E8 to
rank>=9. It does NOT seal the removed own2 automatically. Its own
0/3/4 alternatives are full-star affine E8/cycle/star; own5 retains
the nonpositive3,5,3 path after deletion3. For r>=10 the right D
extension already contains affine E8 and cannot be finite. There
are no other finite right contacts. Thus the common2/3 inventory is
complete with arbitrary spectators, including the exceptional E8 case.

### B.2. The actual pressure packages and their scope

The preceding finite lists reduce distinct1/2 and2/3 to four explicit
pressure calculations. Every spatial vector below is supported on the
same first seven coordinates in B and D, and the actual combinations
use only c and n1,n2,n3. This is why the numerical pressure packages
transfer between B and D; no D-tail statement is imported into B.

* An actual edge2 with cosh>=1+sqrt2 gives h2^2>=sqrt2-1. Distinct6/7
  requires h6+h7>3/2. Orthogonal Cauchy about u2 has projection
  rho=1/sqrt3+sqrt(2/7)<1.114 and residual square
  29/21+4sqrt(2/21)<2.62. Its decreasing endpoint envelope bounds
  h6+h7<1.494416, a contradiction.
* Isolated5@1,Q0 actually also closes4, giving h4^2>=phi-1. For
  w=u6+u7, rho=w.u4<787/500 and perpendicular square<69/50.
  With sqrt(9/14)<401/500 and9/14-(phi-1)<1/40,
  x.w<362087/250000<3/2, contradicting the same distinct6/7 law.
* Isolated5@2,Q0 gives x1>31/50 and S3>131/100. Distinct4/5 demands
  x.w>57/40 for w=u4+u5. Put v=prefix3, k0=1/2+1/sqrt5,
  x0=(31/50,69/200,69/200,0,...), and
  W=x0+(w-k0*v)/6. Its old e1,v,w coefficients remain positive.
  The residual square is<11/10 and its extra demand>9/50. The
  exact Cauchy lower seal at the conservative energy31/50 is

      (31/50+(9/50)/6)^2
         -(9/14)*(31/50+(11/10)/36)=3/700>0.

* For4@2+3@3, the actual white A2 cost is2/3, so Q0 or one A1_3.
  At Q1/4 the actual edge has cosh L3=4+3sqrt2>7. At Q0 the two
  edge floors have equality-vector energy E0=(6sqrt2-4)/7>.6407.
  The positive old-A2 correction makes an ACTUAL future point of
  norm -sqrt2/3 and displacement

      U=((2+sqrt2)/2)e1+((1+sqrt2)/3)*(e2+e3+e4).

  After orthogonal projection onto span(u2,u3), the residual square
  is<1.055 and the pressure gain>.06. Both old dual coefficients
  stay positive on adding this residual with coefficient1/16. Thus

      (.6407+.06/16)^2
        -(9/14)*(.6407+1.055/256)
          =8825113/11200000000>0.

The last package uses an actual timelike point, not just two edge
floors. Replacing it by an edge-covering claim would be false.

The remaining isolated7@1 has a different actual point. Set
t=cos(2pi/7), g=cos(pi/7), and

    c=sqrt(t)*p-g*sqrt2*e1,
    xi=c+n1/g+(1-t)*n2/g,
    V=t(e1+e2)+(1-t)e3.

All coefficients are positive, xi.c=xi.n1=0,
xi.n2=(1/2-t)/g<0, and

    xi^2=-(1-t)(2t-1)/(1+t),
    xi=sqrt(t)*p-V/sqrt(1+t).

It is future and feasible against all actual walls. The cubic bracket
.6234<t<.6235 gives
x.V>=sqrt(t(1+t))-sqrt((1-t)(2t-1)) and x.(V/||V||)>.73.
For w=u6+u7 its projection onto this unit direction has1.3<rho<1.35
and residual square<11/5. The Cauchy envelope decreases for pressure
>=.73, giving

    x.w<1.35*.73+sqrt((11/5)*(9/14-.73^2))<591/400<3/2.

The final radical margin is277/1750000. This closes the last1/2
role. Every auxiliary point is an actual positive combination, so
arbitrary unknown walls and spectators remain satisfied.

With distinct1/2,2/3,6/7 established, the manuscript's fixed three-pair
dual proves the actual9/14 floor. It is not necessary to duplicate
that already displayed Gram sum here.

## C. The large caps: clearing, finite continuation and positive duals

The cap proofs have three logically separate methods. Keeping them
separate prevents a whole irreducible vertex argument from being
misapplied to a component with spectators.

### C.1. Clear spectator pressure BEFORE transporting a component proof

For an isolated own contact gamma at i, positively clear the whole
foreign forest. Write d0=1-Q>0, u=d0/gamma^2 and
A^2=gamma^2(D_i-u). If0<u<1 the explicit positive combination

    xi=c+v_foreign+(d0/gamma)*n_i

has products0 with c and every foreign old node, negative product
with n_i, and norm -gamma^2*u*(1-u). It is future and feasible
against every actual wall. Ifu=1 the same positive vector is null,
forbidden compactly. Ifu>1, the actual positive set consisting of
c,n_i and the foreign old nodes gives the finite-projection pressure.
Together these give the PURE component inequalities

    h_i>sqrt(D_i-u)-sqrt(u(1-u))     (u<1),
    h_i>sqrt(D_i-u)-sqrt(u-1)        (u>1).

In A the sharper finite form is
h_i>sqrt(D_i-1)-sqrt(u-1), using x.n_i>0. These are not projected
facet normals asserted native. With Q>=1/4 and gamma>=phi/2, they
make every foreign isolated central cap strong at the .61 loss.
The u<1 scalar maximum is elementary; the u>1 bound uses
u<=3/phi^2 and the displayed lower range of D_i.

For a zero-own far3 at endpoint1, B/D have A^2=Q-1/2<1/(2i),
while a positive foreign contact has derivative square>=1/4.
It would therefore close an orthogonal component, impossible. In A,
zero-own far3 at1 has height<1/(2m), and next-to-end3 has
Q<2/m<=2/5 and A^2=Q-2/N<Q. Thus both are excluded by(A2);
an endpoint4 has Q<1/m and negative height. These repairs are
essential: positive full-star rank alone no longer removes a
zero-own wall when foreign components exist.

For a foreign far3 helper, the residual drops by Q>=1/4. The exact
edge ratios are uniformly above the weak target on the relevant
finite index range5..29. For A, after writing the white path parameter
m>=5 and f=101/125, the bound is

    (2fm+1)/sqrt(m(m+2))>=227/(25sqrt35)>761/500.

The derivative sign is(2f-1)m-1>0. The native near4/simple3 foreign
role has ratio squared>121/35>(761/500)^2. Thus every weak A cap
has Q0. The analogous B/D foreign-far bound uses
cosh L_i=(2i*gamma+1)/sqrt(i+1-2iQ), and retains Q throughout.

### C.2. Finite continuation lists, including exceptional and zero-own roles

At a central release the left old A path has complete one-node finite
extensions endpoint3,endpoint4,next-to-end3, or an E6/E7/E8 position3.
Their residuals, with path parameter i, are

    (i+1)/(2i), 1/i, 2/i, and1/4,1/7,1/16.

Two contacts to the same retained piece make a cycle. A full bridge
joining both old pieces has large rank; finite classification forces
the corresponding A/B/D simple endpoint bridge. Its new node is
internal, so a foreign extension would create the wrong branch or
an internal4 in a finite diagram. This is a proved local seal, not
a general statement that E7, D7 or another intermediate piece has
no extension.

The right finite possibilities used by the caps are exactly:

| Retained right piece | Extra finite one-node attachments |
| --- | --- |
| B_s,s>=4 | simple3 at long end |
| B3 | long3, or short3 giving F4 |
| B2 | simple3 at either end |
| D_s,s>=5 | long-leaf3; fork3 to E_(s+1) only s5,6,7 |
| D4 | simple3 at ANY leaf |
| D3=A3 | middle3, end3, end4, end5 |

This table is derived by spherical finite classification or the leaf
inverse-Gram costs; multiple contacts make a cycle. The D fork E6/E7/E8
roles are RETAINED with every foreign extension allowed by the full
finite replacement. They are not sealed merely because their first
attachment is exceptional.

In B the short-B3 extra3 seals to F4; own0 is affine and own>=3
is strong. The short-B2 extra3 has own0 affine, own4 affine star,
own6 affine3,6 and own>=5 strong. Its sole additional weak role is

    Q0, v=(omega_(r-2)+omega_r)/2, own3.

If it joins a left endpoint too, the full white B_r bridge is strong
even at own0 and has no foreign extension. The weak own3 role has
an ACTUAL strong r-1 mate, from its positive A_r replacement, with

    cosh L_(r-1)=[(1+sqrt2)*(r-2)+1]/sqrt(r+1).

It is not legitimate to remove this role by an F4-sealing assertion.

For the D4 short layer i=r-4, fork3 has
Delta=1/2-Q, sigma=gamma+1/4. Own0/3/4/6 are exact affine
principals; Q>0 and own>=7 are strong. At Q0 only own5/fork3
needs the averaged point below. The large endpoint bridge is D_r
and strong (own0 gives an affine cycle).

For D3 at i=r-3, middle3 is strong or affine. End4 has
Delta=1/4-Q and is strong once own0's affine fork is removed.
End5 is sealed H4 and is strong EVEN AT own0. End3 has
Delta=5/8-Q and sigma=gamma+1/4. Its left endpoint bridge is A_r,
sealed against further foreign contact and strong; own0 is a cycle.
Without that bridge, own0/3/4/6 are affine, Q>0 is strong, and
gamma>=1 is strong. The only extra roles are own5 or finite>=7,
gamma<1, with fork3, Q0 and no left contact. This is why the D cap
needs two distinct actual averaged feet, rather than treating D3
as an ordinary long tail.

### C.3. Actual mixed replacements and the endpoint certificate

In B/D all remaining ordinary weak indices have exactly

    v=gamma*omega_i+omega1/2,
    A^2=2i*gamma^2+2gamma-1/2,
    Q0, 101/125<gamma<6/5.

For adjacent weak i,j, weighted AM--GM of their heights bounds their
wall product strictly above -1/2; the native gap makes it0. After
retaining all old walls except i,j, the two-cap Schur matrix is
I-t*ones(2,2), t=(i-1)/(2i), positive definite with positive inverse.
It therefore defines an ACTUAL mixed endpoint. Its positive pressure
exceeds what two simultaneous .61 losses permit: the total loss is
<7/16+1/(4sqrt10)<.61. Thus ordinary weak indices are nonadjacent.
All inverse products are C_ij=2min(i,j), so this argument takes place
in the actual component and appends the full spectator Gram unchanged.

At own5/fork3 in D4/D3 put gamma=phi/2. Two ACTUAL H4 projections,
one on either side of i, have positive corrections
gamma*(3n_i/2+n_side+n_second/2), residual delta5=1-3gamma^2/2.
The right second node at the D3 layer is the OTHER, uncontacted fork.
Their average Y has

    x.Y>A-sqrt(delta5), Y^2=A^2+1-15gamma^2/8<A^2.

At i>=19 this gives normalized demand>.97 and correlation with u5
<.55. Strong5 demands>.8; their sum exceeds1.77 whereas its squared
norm is<3.1<1.77^2. For finite own>=7 at the D3 layer use the two
negative old-A2 cores with corrections
gamma*(4n_i/3+2n_side/3). Their average satisfies

    x.Y>=A-sqrt(d), d=4gamma^2/3-1,
    Y^2=A^2+1-14gamma^2/9<A^2.

At i>=21 its correlation with u5 is<1/2 and d/A^2<=1/195,
so its normalized demand>.928; strong5 gives>.807. Their sum1.735
has square3.010225>3, a contradiction. These are exactly the accepted
D24-and-higher tail bounds. NO pending D23 whole-vertex extension
is used in this appendix.

Index5 is strong before either tail argument: if weak its actual
finite/core/dotted feet force x1>.5965, while6 and one of16/17
are strong, giving a positive-prefix interpolation of energy>1.
This makes the tail assembly noncircular. If21 is strong, the
concave prefix reference has E21>1.009 and at most eight isolated
holes cost<8/3200. If21 is far-weak,20 is strong and its actual feet
give x1>.424, yielding

    E20+(5/4)*(.424-.361)^2-7/3200
       >.99752+.00496125-.0021875=1.00029375.

For B the one additional terminal weak role is handled by its strong
mate and the accepted isolated/two-hole terminal loss bounds. This
closes B23..29 and D24..29 with arbitrary spectators.

For A21..29 the same METHODS apply but the exact endpoint certificate
is different. Its central inverse is
C_ij=2min(i,j)(N-max(i,j))/N, h0=hN=0, and
||x_A||^2=(1/2)sum(h_(i+1)-h_i)^2. After C.1 clears all foreign
pressure, weak central caps are near4/simple3 or literal far-end3
helpers with own5,finite>=7 or dotted. The near4 mate has normalized
height>sqrt3/2 and is excluded against an opposite strong anchor.
For literal packets, same-helper native products force0 and
opposite-helper products force -1/2; positive actual mixed replacements
then forbid adjacent weak packets. These product statements require
the complete two-helper AM--GM estimates, not an arbitrary same-helper
packing rule.

The A21 reference and its weak5 duals lift to N22..30 by inserting a
constant potential plateau. Its potential norm gains only zero
differences, so the lift is ISOMETRIC; the positive central coefficient
splits between the two plateau ends. The exact opposite-foot projections
remain on the unchanged end windows. Strong central floors, removal of
at most six holes and the endpoint certificate give energy>4013/4000.
An increase of the scalar floor alone would not prove this transfer;
the full positive-dual potential lift is required. This completes the
grouped account of the three cap methods, subject to the transcription
boundary below.

## D. What is ready to insert, and what remains companion-dependent

The appendix has deliberately separated proof from audit provenance.
The ordinary rates and short-debit analysis in A.1--A.7 are expanded
here as mathematical arguments; their source versions are O1--O8 in
the existing manuscript map. No new hash table or new hypothesis is
needed in the mathematical body. A finite-type or Lanner classification
is an explicit external input; every component energy estimate remains
an internal proof.

The following portions of B/C are coherent, prose-ready outlines with
their principal formulas, but are NOT yet a complete replacement for
every companion scalar derivation:

1. The sharp scalar envelopes in B.2 for the edge2 threshold and the
   two isolated5 pressures need their full monotonicity and positive
   rational seals inserted from the low-index pressure source.
   The native7 cubic/normalization and decreasing-Cauchy seals must
   accompany its displayed actual point; the point alone is insufficient.
2. The complete classical bulk left-position inequalities, exact
   isolated-cap feet and foreign-far target comparisons are summarized
   in C.1--C.2. Before removing the current companion boundary, include
   their full finite residue/derivative tables from C1--C3 and P1--P2.
   Do not cite only the output weak profile as a classification theorem.
3. The B terminal-role strong-mate half-tangent inequality and the
   one terminal two-hole loss still need their full rational estimates
   in a standalone version. The D two averaged-point products are
   displayed, but their exact elementary radical seals should travel
   with them, rather than being replaced by decimal evaluations.
4. The A21 plateau construction requires its explicit potential
   coefficients and the complete same/opposite-helper determinants,
   near4 anchor choices and endpoint-hole argument. The last paragraph
   of C.3 is a structural proof map, not those missing algebraic pages.

These derivations are present in the internally reviewed companion
packet; the work still pending is their complete integration into a
standalone exposition and review of that translation. Until then retain
the primary abstract/scope qualification. No source hash or review PASS
turns these finite inventories into published assumptions, and no new
computation is needed merely to transcribe their existing proofs.

The existing global23 composition, accepted A20>11/14 and accepted
D23-plus-isolated-I argument are the ONLY endpoint inputs here. The
held new whole-D23, B22/D22-plus-I and uniform A4/7 candidates are
excluded from this staging note. Root owns integration into the ONE
compact manuscript, mathematical review and any later bounded build.
