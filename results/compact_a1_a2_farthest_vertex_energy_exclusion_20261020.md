# A1/A2 forests cannot be farthest finite vertices in dimension at least14

20 October 2026. Complete paper argument; no execution, enumeration,
optimizer, build or formal verification. Review is pending. The argument
uses no arithmetic hypothesis.

**Theorem.** Let P be any compact hyperbolic Coxeter n-polytope, n>=14,
and choose any interior future unit point z. A vertex p maximizing
E=-z.p over the actual vertices cannot have whole diagram A1^a+A2^r.
Here a+2r=n. Other vertices of P are unrestricted.

This is a POINTWISE restriction at the one farthest vertex, not a new
dimension ceiling for all compact Coxeter polytopes. Its intended use is
to remove this entire finite-diagram family from a future uniform bound.

The farthest-vertex mechanism and shadow inequality were developed with
batch_guard; the present proof supplies the full A1/A2 inventories and
energy bookkeeping. It combines and strengthens the two reviewed support
notes `compact_right_corner_short_edge_capacity_20261020.md` and
`compact_a2_corner_short_edge_capacity_20261020.md`. All local interfaces
needed for this theorem are derived below.

## 1. The pointwise energy and actual-wall shadow inequality

Use outward unit normals n_i at p, with finite Gram B and dual vectors
omega_i. Put d_i=(B^-1)_ii and u_i=omega_i/sqrt(d_i). Write

    z=E p+z_T,       x=-z_T/E,       ||x||^2=1-E^-2<1.

If the outgoing i-edge has length L_i and opposite endpoint p_i, then
p_i=cosh(L_i)p-sinh(L_i)u_i. The inequality -z.p_i<=E gives

    x.u_i >= tau_i=tanh(L_i/2).                         (1)

The orthogonal A1/A2 summands of p^perp have energies summing to ||x||^2.
An A1 energy is at least tau_i^2. In an A2 block u_1.u_2=1/2; thus its
energy is at least either individual tau_i^2, and also at least

    (tau_1+tau_2)^2/3.                                 (2)

For reference the exact constrained minimum is
4(tau_1^2+tau_2^2-tau_1*tau_2)/3 when the ratio is in [1/2,2],
and otherwise the square of the larger constraint. Only (2) and the
individual constraints are required here.

For an actual outside wall c, write contacts beta_i=-c.n_i>=0 and
height A=-p.c>0. With t=B^-1 beta, its spatial square is

    beta^T B^-1 beta=1+A^2.

If c itself closes edge i, eliminating the white normals gives

    delta=1-beta_W^T B_W^-1 beta_W>0,
    A^2=delta*sinh^2 L_i,       delta<=1.                (3)

There is also a useful bound when a DIFFERENT wall d closes that edge.
Let its contacts be gamma, height D, and q=B^-1 gamma. The inverse of a
finite Coxeter Gram matrix is entrywise nonnegative. The Schur identity

    (B^-1)_jk-(B^-1)_ji(B^-1)_ik/d_i
        =(B_(I\{i})^-1)_jk >=0       (j,k!=i)

therefore proves, for beta,gamma>=0,

    beta^T B^-1 gamma >= t_i*q_i/d_i.

All-wall acuteness for the DISTINCT actual walls gives
A*D>=beta^T B^-1 gamma. The actual endpoint equation for d gives
tanh L_i=D*sqrt(d_i)/q_i. Consequently

    tanh L_i >= t_i/(A*sqrt(d_i))       if d!=c.         (4)

The distinct-wall hypothesis will always be justified explicitly. No
claim that reflections of a facet form a discrete group is involved.

We also use actual finite endpoint recognition: if c plus the n-1
retained old normals is positive definite, projection of p to their
common orthogonal line gives the actual adjacent finite vertex. Indeed
its correction coefficients are nonnegative by the same inverse-Gram
property, so all other actual halfspaces remain satisfied by acuteness.
The omitted old wall is strict in each application below. Thus a formal
positive-definite corner is not silently substituted for an actual one.

## 2. Numerical constants are rational comparisons

Set

    L*=asinh(1/sqrt2),       tau*=sqrt3-sqrt2,
    sinh^2 L0=31/50,         tau0=tanh(L0/2),
    e=319/2000.

Then

    tau*^2>1/10,       tau0>173/500,
    (4/3)tau0^2>e.                                    (5)

These follow by squaring radicals and using
sinh^2 L=4tau^2/(1-tau^2)^2. The last comparison is
4*173^2/(3*500^2)>319/2000.

Put phi=(1+sqrt5)/2 and q=cos(2pi/7). The elementary bounds

    623/1000<q<5/8,       4sqrt2/9>31/50                (6)

follow from the increasing polynomial8q^3+4q^2-4q-1 on [1/2,1]
and rational squarings. In particular cos^2(pi/7)>1623/2000.
We will use these exact half-tanh bounds for the isolated7 continuation
at an A2 block:

    tau_contact>7/50,       tau_mate>3/10.              (7)

Its sinh-squares are (2q-1)/3>41/500 and
2(2q-1)/(3(1-q))>164/377. Substitution proves (7); for the second
comparison164*8281-3600*377=884>0.

Every actual nonzero contact is a native cos(pi/m), m>=3, or exceeds1.
Value1 and positive affine principal subdiagrams are forbidden by
compactness: their positive null combinations are future feasible ideal
points, using all-wall acuteness. In particular the affine3--6 path,
4--4 path,3/3/4 star and four3-leaf star cannot occur. These elementary
null diagrams will exclude specific cases below.

## 3. Complete A1-short inventory and its energy bound

Call an A1 edge short when L<L*. By (3) any closing wall has A^2<1/2.
Two DISTINCT such walls have product

    c.d=beta^T B^-1 gamma-A D>-1/2,

so are orthogonal by the Coxeter gap. An A2 block touched at both nodes
would have white cost at least1 when an A1 edge is released, impossible
for a finite endpoint. A single contact beta in an A2 block contributes
(4/3)beta^2, at least1/3.

First suppose c has only A1 contacts. Write their squares s_j. The
conditions are sum s_j>1 and, at a released contact i,

    2s_i+3sum_(j!=i)s_j<3.

There are at most three support nodes. Four force s_i<3/8, hence s_i=1/4;
the other three must also be1/4, giving height zero. For three support
nodes the released square is below3/4. The possibilities1/4 and1/2
would require the other two squares respectively in (3/4,5/6) and
(1/2,2/3), neither a sum of two nonzero native squares. The remaining
square cos^2(pi/5) requires both other contacts3. With two support nodes,
the smaller square is1/4 or1/2: any larger one makes five times that
square exceed3. The resulting complete list is

* isolated dotted contact, with height squared in (0,1/2);
* (3,m), finite m>=7, with height squared in (0,1/4);
* (3,dotted), with height squared in (1/4,3/8);
* (4,5), or (3,3,5), with height squared H=(sqrt5-1)/8.

For (3,m) only7 can close two short directions: the second direction
requires cos^2(pi/m)<5/6, while label8 already exceeds that. The (4,5)
wall has two short directions; (3,3,5) only its5 direction. An isolated
dotted wall and a (4,5) wall cannot coexist with any other wall in this
A1-only list: any shared coordinate already contributes more than the
product of their heights. For (3,3,5), orthogonality could only use one
3--3 overlap, forcing the other squared height1/(16H)>3/8, whereas
every listed wall with a3 contact has height squared below3/8.

Only finite (3,m) and (3,dotted) can coexist. Their common coordinate
must be the light3 contact on both, giving height product1/4. There is
at most one of each type. If the finite wall closes two short directions
(label7), its height squared is below1/12; this is incompatible with
the dotted partner's squared height below3/8. Thus there are at most
two short A1 directions in the A1-only subfamily.

Now allow A2 contacts, and let s be the released A1 contact square and
Q the full white cost. Shortness says

    1-s<Q<1-(2/3)s.

For released labels3,4,5,6 the intervals are respectively
(3/4,5/6),(1/2,2/3),(1-s_5,1-2s_5/3),(1/4,1/2).
The native cost gaps1/4,1/3,1/2,2/3 leave, with at least one A2 contact,
only Q=7/12 for label4 or Q=1/3 for label6. The first is the affine
3/3/4 star and the second the affine3--6 path. For label5, the next A2
cost after1/3 is2/3, and adding even an A1_3 costs7/12, above its upper
bound;1/3 is below its lower bound. For label3 the only intervening
candidate before5/6 would have to refine native gaps: one A2_3 plus
an A1_4, or two A1_3, is already5/6, while the smaller combinations
are at most2/3. Thus no omitted case lies in its open interval.

For released finite m>=7 the upper bound is below1/2. The only possible
A2 contribution is a SINGLE3 contact, with no other white contact.
A released dotted contact gives an upper bound below1/3 and is
impossible here. Therefore the ONLY mixed A1-short profile is

    M_m=(A1_m,A2_3),       finite m>=7,
    A^2=cos^2(pi/m)-2/3<1/3.                            (8)

Such a wall has exactly one A1 direction. It cannot coexist with an
A1-only short wall: a shared A1 contributes at least sqrt3/4, greater
than their height product (<sqrt(1/6)); disjoint support would give a
product strictly between-1/2 and0. Two M walls must have different A1
coordinates and use OPPOSITE nodes of the SAME A2 block: same-node
overlap contributes1/3, and different blocks give no spatial overlap.
Hence at most two M walls coexist. No fine native-pair equations are
needed to obtain this count.

This proves at most two short A1 directions in total. If there are two,
their squared half-tanh sum exceeds2/25. For the single (3,7) wall the
bounds are tau_heavy>7/50 and tau_light>13/50, whose squared sum exceeds
2/25: its light-edge sinh-square is
(2q-1)/(2(1-q))>123/377, which implies the latter bound. For (4,5)
both half-tanhs exceed1/4. A finite/dotted two-wall pair
has squared heights respectively>1/6 and>1/4, giving half-tanhs>11/50
and>13/50. Finally each M direction has
sinh^2 L>(3/2)(1623/2000-2/3)=869/4000, giving tau>11/50.
All comparisons are direct substitutions in the identity preceding(5).

Using (5) on the remaining A1 directions, in EVERY case their total
energy therefore satisfies

    E_A1 > a/10-3/25.                                 (9)

If there is no short A1 direction, the stronger E_A1>=a/10 holds,
strictly when a>0.

## 4. Complete A2-short inventory at the longer threshold

Call an A2 direction short here when L<L0. Its closing wall has A^2<31/50.
It cannot touch both nodes of an A2 block: two3 contacts already form
an affine triangle. The next possibility(3,4) gives
A^2>=sqrt2/3 and delta<=3/4 at its shorter release, hence
sinh^2 L>=4sqrt2/9>31/50. Releases outside that double-contact block
leave white cost at least1 and cannot be finite.

Thus each touched A2 contributes q_j=(4/3)beta_j^2. Consider the A2
block containing the released node. If c touches its mate instead,
releasing the contacted node is strictly shorter, so it suffices first
to classify the contacted release. Let B0=beta^2 in this block and Q
the complete off-block white cost. Then

    A^2=(4/3)B0+Q-1,
    max(0,1-(4/3)B0)<Q<1-(200/243)B0                  (10)

when there is off-block contact. Each off-block A1 contribution is a
native square, each off-block A2 single contribution is4/3 times one.
A double A2 contact cannot be in this positive-definite white support.

For beta label3 the window is(2/3,193/243). Its only native sums are
3/4: one A1_6, A1_3+A1_4, or three A1_3. They give respectively affine
3--6,3/3/4 and four3-leaf stars. For label4 the window is(1/3,143/243).
The possibilities are1/2 (A1_4 or two A1_3), or7/12 (A1_3+A2_3).
All contain affine4--4 or3/3/4. These lists are complete because any
further contact costs at least1/4, and the next single costs are1/2
and2/3. Label5 leaves just Q=1/4 or1/3. Label6 is forbidden by its
old A2 mate and the affine3--6 path. For labels>=7, (10) leaves only
Q=1/4; Q=1/3 is already too large by q>31/50 in(6). This is short
exactly when B0<729/800, so only labels7,8,9,10 can occur. Indeed
cos^2(pi/10)=(5+sqrt5)/8<729/800, while
cos^2(pi/11)>1-(2/7)^2=45/49>729/800, using pi<22/7.

The COMPLETE list of walls closing an A2-short direction is therefore

    I: one isolated A2 contact, finite m>=7 or dotted;
    F: (A2_3,A2_5);
    N_m: (A1_3,A2_m),       m=5,7,8,9,10.              (11)

An untouched A2 block cannot supply a finite endpoint, since its white
cost would be the full spatial square>1. The classification therefore
also covers all releases with zero contact at the released node.

## 5. Every A2 block pays at least e, except one isolated7 block

If both of a block's directions are at least L0, (2) and(5) give energy
greater than e. For a short direction use(11).

For finite isolated I, c plus the retained old normals is positive
definite at the mate release, so the mate endpoint is ACTUAL by the
projection argument in Section1. With B0=cos^2(pi/m),

    sinh^2 L_mate=[(4/3)B0-1]/(1-B0).

For m>=8 this is at least2sqrt2/3, giving tau_mate>2/5 and block energy
>4/25>e. For m7, (7) gives block energy>9/100. Call this the exceptional
isolated7 block.

For dotted I the mate white cost B0>1 is not positive definite, so its
actual opposite is DISTINCT from c. Apply(4) with d_i=4/3 and
t_i=2beta/3. Since A^2<31/50 implies B0<243/200,

    tanh^2 L_mate >= B0/(4B0-3)>81/124>16/25.

Hence tau_mate^2>1/4>e.

For F its5-mate corner is positive definite:
delta=2/3-phi^2/4>0. The ACTUAL mate edge has cosh L=phi^3,
so its squared half-tanh is1/phi>e. For N5 its A2-mate corner is also
positive definite, with delta=(3-phi^2)/4>0. That ACTUAL mate edge has
cosh L=phi^2/sqrt3>3/2, hence squared half-tanh>1/5>e.

For N7..10, set B0=cos^2(pi/m), H=(4/3)B0-3/4. Its A2-mate white cost
B0+1/4>1, so its actual opposite is DISTINCT. Equation(4) yields

    tanh^2 L_mate >= B0/(3H)=B0/(4B0-9/4)>11/17>16/25,

where B0<=cos^2(pi/10)<11/12. Again its block energy exceeds1/4>e.
Thus every case in(11), apart from isolated7, pays more than e.

There is at most ONE isolated7 wall in the whole forest. Two in different
A2 blocks have product minus their common squared height, strictly
between-1/2 and0. Two in the same block have positive spatial product
already greater than1/2, while the height product is below1/12. Both
contradict the actual Coxeter Gram alternatives. Such an isolated wall
cannot close directions in any other block. Consequently

    E_A2>r*e-7/100,                                   (12)

because e-9/100=139/2000<7/100. If no isolated7 occurs, the stronger
E_A2>=r*e holds, strictly when r>0.

Moreover an isolated7 block forbids ALL short A1 directions. An A1-only
short wall has disjoint spatial support and height squared<1/2, so its
product with the isolated7 wall is strictly between-1/2 and0. An M wall
in(8) has height squared<1/3. If it uses another A2 block the same gap
applies; if it uses the isolated7 block, even opposite-node spatial
overlap is beta7/3>sqrt3/6, whereas the height product is below1/6.
The product is positive. Thus (9) improves to E_A1>=a/10 whenever the
exceptional isolated7 block occurs.

## 6. Contradiction to the finite-point energy

If isolated7 occurs, add(12) and the improved A1 bound:

    ||x||^2 > a/10+(319/2000)r-7/100
             >= (319/4000)n-7/100 >1       for n>=14.

If it does not occur, add(9) and the stronger A2 bound. For n>=15,

    ||x||^2 > a/10+(319/2000)r-3/25
             >= (319/4000)n-3/25 >1.

At n14 the same inequality suffices if a>=2, since r<=6 and its right
side is at least1/5+6*(319/2000)-3/25=1037/1000>1.
If a=0, use E_A1=0 instead of(9): E_A2>7*(319/2000)>1.
These cases exhaust a+2r=14. All contradict ||x||^2<1 in(1), proving
the theorem.

No assumption was made about the diagrams at any other vertex, a shared
global choice of white weights, arithmetic fields, or a classification
of finite-volume high-dimensional chambers. Standard all-wall acuteness,
compact simplicity and the actual native Coxeter-angle alternatives are
the only geometric inputs beyond the explicit calculations.
