# B5--B8 component energy with arbitrary finite spectators

20 October 2026. Complete paper candidate, awaiting root and outside
review. No scientific execution, optimizer, Lean, TeX, literal-host
classification, or accepted-source edit. Derived by batch_guard in
discussion with root; the terminal bridge inventory is written here
explicitly rather than transferred from a larger B diagram.

**Theorem.** At a vertex farthest from a strict interior point of a
compact hyperbolic Coxeter polytope, the projection onto any whole
finite B_r component satisfies

    r=5 or6:  ||x_B||^2>3/7;
    r=7:      ||x_B||^2>1/2;
    r=8:      ||x_B||^2>4/7.

In particular each B5--B8 component pays strictly more than rank/14.
All orthogonal finite spectators are arbitrary. No ambient dimension
ceiling is inferred in this source.

## 1. Actual geometry and one-pair forcing

Use unit roots n_i=(e_i-e_(i+1))/sqrt2 for i<r and n_r=e_r.
The spatial dual directions are

    u_i=prefix_i/sqrt(i) for EVERY1<=i<=r,
    omega_i=sqrt2*prefix_i (i<r), omega_r=prefix_r.

Let x_B be the component of the common farthest vector x, ||x||<1.
Its actual finite edge half-tangents satisfy
x_B.u_i>=tau_i=tanh(L_i/2). For DISTINCT actual opposite walls,
reciprocal shadow and all-wall acuteness give

    F(tau_i)*F(tau_j)>=k=u_i.u_j,
    F(t)=2t/(1+t^2).                                      (1)

If ||x_B||^2<=e, Cauchy gives
tau_i+tau_j<=sqrt(2e(1+k)). Since F is increasing and log F is
concave on(0,1), the left side of(1) is at most

    4z/(1+z)^2, z=e(1+k)/2.                              (2)

The following exact comparisons therefore FORCE a common opposite:

* r5/6, e3/7, pair2/3: k=sqrt(2/3). The cleared difference
  k(17+3k)^2-168(1+k)=127k-100 is positive, since
  2*127^2>3*100^2.
* r7, e1/2, pair6/7: k=sqrt(6/7). The cleared difference
  k(5+k)^2-16(1+k)=(69k-52)/7 is positive.
* r8, e4/7, pair6/7: the cleared difference
  k(9+2k)^2-56(1+k)=(199k-176)/7 is positive.

The last two use k>9/10. These are full-tangent products derived
from HALF-tangent energy, with the conversion retained explicitly.

Suppose henceforth that c is this common actual opposite. Write

    c=A*p-v_B-vS, A>0,
    beta_i=-c.n_i>=0, Q=||vS||^2>=0.

Here Q contains the FULL inverse-Gram costs of all foreign contacts.
If Q>0 then Q>=1/4, since every nonzero finite native contact is at
least1/2 and a finite unit-Gram inverse has diagonal at least1.
Both replacements deleting the two closed directions are ACTUAL
positive definite Grams. Thus every contact retained by either one
is finite native, and c contacts at most one node of each connected
retained old component. In particular multiple contacts with a path
would make a nonpositive cycle and are not permitted in that finite
replacement.

For an actual direction i, let a_i=v_B.u_i. If the retained white
residual is delta>0, then

    cosh L_i=a_i/sqrt(delta),
    tau_i^2=(cosh L_i-1)/(cosh L_i+1).                     (3)

Subtracting Q reduces delta and only lengthens that edge. In the
occasional extra direction below, a_i^2>A^2 forces c itself to be
its actual opposite by reciprocal shadow; positive acute projection
establishes actualness, rather than assuming a formal completion.

The generic geometric inputs are the unchanged sources
`compact_farthest_vertex_energy_obstruction_20261020.md`, SHA256
`50c64e5a521d3e4e0f277a8ece02d3a29aa4ae8b71485cf3088211e587f45bff`;
`compact_actual_wall_reciprocal_shadow_20261020.md`, SHA256
`c913ec3bd8e0f26a8c679a1ede2b98b251ad3f36efeae60486213ddfa3e8badc`;
actual-point source SHA256
`3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`;
and orthogonal-component opposite obstruction SHA256
`4216d06c2e03fd336898b8dac289b9969fa096c8a1d08bd67ddeb8ab5ec2f23d`.
We use finite Coxeter classification and the external classical fact
that a compact Lanner diagram has at most5 nodes. A compact chamber
has no affine principal diagram: its positive null combination would
be an actual feasible ideal point. All paths below are exact selected
principals; contacts outside their displayed nodes cannot remove them.

## 2. B5 and B6: complete common2/3 inventory

Suppose ||x_B||^2<=3/7, so c closes2 and3. Deleting2 retains old1
as A1 and a right B_(r-2). For r6 this is B4: its only finite
one-node extension has simple3 at its long endpoint3. For r5 this
is B3: it additionally permits simple3 at its short endpoint5,
making F4. No other right contact is finite.

First consider no short-F4 contact. If c contacts3, that contact is
3. Deleting3 retains left A2, so c has only one left contact, at1
or2, of label3/4/5. A zero contact at2 makes the fork with leaves
(old2,c) and right arm ending in4 an exact affine diagram, even if
c also contacts old1. Thus c contacts2 and not1. Own3 is the affine
triangle c,2,3. For own gamma>=cos(pi/4), the left white cost is
4gamma^2/3; hence

    delta<=1-4gamma^2/3<=1/3,
    a_3=(4gamma+3)/sqrt6,
    cosh L3>=2+3/sqrt2>3.                                (4)

All Q are included by the first inequality. Equation(3) gives
tau3^2>1/2>3/7, a contradiction.

With no right contact, c has no B derivative unless it contacts old1
or2. The COMPLETE possibilities, by the retained A2, are one contact
3/4/5 at either node. A simple3 at2 gives the affine fork ending in4;
4 at1 gives the double4 affine path;4 at2 gives the affine(3,3,4)
star on old1,old2,old3,c. These are impossible independently of Q.

For3 at1, A^2=Q-1/2>0, while actual deletion3 has left cost1/3,
so Q<2/3 and A^2<1/6. Every contacted spectator node has squared
derivative at least1/4>A^2, and therefore is also closed by c.
This contradicts the farthest orthogonal-component opposite theorem.
If there is no spectator contact the height is already negative.

For5 at1 or2, put phi=(1+sqrt5)/2 and gamma=phi/2. Actual deletion3
has left cost phi^2/3>3/4. Thus Q<1/4 and consequently Q0.

For5 at1, A^2=(phi-1)/2. At direction4,

    a_4^2=phi^2/8,
    a_4^2-A^2=(5-3phi)/8>0,
    cosh L4=phi^3>3.                                    (5)

The positive difference supplies the ACTUAL extra closure4 via shadow.
Hence tau4^2>1/2, a contradiction. For5 at2 the already actual
direction3 has

    cosh^2 L3=2phi^2/(2-phi)>9,                          (6)

where the last comparison is11phi>16. This is another contradiction.
All no-right profiles and foreign costs have now been covered.

Finally retain the r5 short-F4 possibility: c contacts5 by3 while2
is deleted. This F4 seals all contacts with old1 and every spectator,
since no larger finite connected diagram contains F4. Thus only the
removed2 can be additionally contacted. Its zero contact gives the
affine F4 path old2,old3,old4,old5,c. In the actual deletion3,
the connected white path old1,old2,c,old5,old4 has rank5 and its
terminal4. Finite classification forces c--old2 to be3: label4
would give the affine4,3,4 subpath, and a5 or higher cannot occur.
For this sole profile,

    delta=1/4, a_2=1+1/sqrt2,
    cosh L2=2+sqrt2>3.                                  (7)

Again tau2^2>1/2. Thus no common2/3 exists under the proposed upper
energy, contradicting Section1. Both B5 and B6 pay strictly above3/7.

## 3. A no-right-contact lemma for the terminal pair6/7

This paragraph applies to both B7 and B8. If c has no contact with
the retained right part while6 is deleted, then it has no contact
with old7 (or with old8 in B8). Deleting7 retains the left A6, so
the only possible B contacts are one3 at a node j<=6 or one4 at an
endpoint j=1,6. If there is no B contact, c cannot close either B
direction.

For3 at an internal j>=2, select the affine fork with leaves
(old(j-1),c) and the old right path terminating in4. This includes
the affine(3,3,4) star when the long arm is shortest. For4 at1 use
the double4 affine path; for4 at6 use the affine4,3,3 star in B8,
or the affine4,4 path in B7. All these exact principals are forbidden.

Only far3@1 with arbitrary spectators remains. The full star of all
old walls plus c is indefinite because A>0. Any minimal nonpositive
principal is Lanner (affine is forbidden). It contains c and BOTH6,7,
because each actual replacement deleting6 or7 is positive definite.
Within the B part c connects only at1; old spectator components have
no B edges and attach to it only through c. Therefore a connected
Lanner containing7 and c must contain the whole path c,1,...,7,
already8 nodes. This contradicts the Lanner bound5. Spectators do
not supply a shortcut in that graph. The no-right-contact branch is
therefore impossible, without setting Q to zero.

## 4. B7: all terminal common6/7 profiles have cosh L6>3

Suppose ||x_B||^2<=1/2. By Section1 one actual c closes6 and7.
Deleting6 retains the left A5 and the short old7 as A1. Its contact
beta=-c.n7, if nonzero, is any FINITE native cos(pi/m), m>=3.
Section3 covers beta0. Deleting7 retains A6, so c has at most one
left contact, either3 or an endpoint4 (or the permitted internal3
exceptional positions); label>=5 cannot extend A6 finitely. Write
gamma=-c.n6.

If beta>=cos(pi/5), a left contact j<6 would connect A5,c,old7
into a finite rank7 diagram containing a5 or higher bond, impossible.
Thus no such contact occurs. The actual direction6 has

    delta=1-beta^2-Q,
    a_6=sqrt12*(gamma+beta/sqrt2),
    cosh L6>=sqrt6*cot(pi/m)>=sqrt6*cot(pi/5)>3.           (8)

For the last seal, cot^2(pi/5)=phi^2/(4-phi^2)>3/2 since
5phi^2>12. The lower bound permits gamma0 and every allowed Q.

If beta=cos(pi/4) and gamma0, the selected c--4--old7--4--old6
path is affine. If gamma>0, the retained A6 forces it to be at least
1/2 and forbids every other left contact. Formula(8)'s Schur form
then gives cosh L6>=sqrt24>3. Again all Q only increase it.

It remains beta1/2. If gamma>0, there is no other left contact,
gamma>=1/2, and

    delta<=3/4,
    cosh L6>=sqrt12*(1/2+1/(2sqrt2))/sqrt(3/4)
              =2+sqrt2>3.                              (9)

Let gamma0 instead. If c has no left contact, the path
old4,old5,old6,old7,c is exact affine F4. If its one left contact
is at j1, the A5 endpoint cost forces that contact to be3: a4 would
cost at least5/6, in addition to the1/4 cost at old7. The same
affine F4 principal, omitting old1, remains exact.

For ANY left contact at2<=j<=5, retain its full cost and Q. Its
native contact is at least1/2, so

    delta<=3/4-j(6-j)/12<=1/3,
    a_6>=sqrt12*(j/12+1/(2sqrt2)),
    cosh L6>=1+3/sqrt2>3.                              (10)

The inequality3/sqrt2>2 follows from9>8. Larger left contacts or Q
strengthen this bound; a nonpositive residual is already impossible
at the actual endpoint. Thus(10) retains all finite left bridges,
including the E7 positions, without assuming they are sealed.

Equations(8)--(10) exhaust beta and all zero-own, left and spectator
possibilities. They give tau6^2>1/2 by(3), contrary to the proposed
component energy. Hence B7 pays strictly above1/2.

## 5. B8: the short B2 bridges are retained explicitly

Suppose ||x_B||^2<=4/7. Again c closes6 and7. Deleting6 retains
right B2=(old7--4--old8). Its COMPLETE finite one-node extensions
are simple3 at either endpoint, making B3; a4 attachment is affine
and a larger attachment is nonpositive. Contacting both endpoints
would give a nonpositive triangle. The zero-contact branch is
Section3.

For simple3 at LONG endpoint7, let gamma=-c.n6. If gamma0, old7
with the two simple leaves(old6,c) and its4 leaf old8 is an exact
affine star; any left contacts outside that principal are irrelevant.
Thus gamma>0. The retained A6 while7 is deleted allows no other
left contact and forces gamma3 or4. Own3 makes the affine triangle
c,old6,old7; the only possible case has own4. Its actual direction7
has the standard prefix value

    cosh L7>=6+7/sqrt2>7,                              (11)

with equality at Q0 and larger value for Q>0. This alone violates
the energy4/7 upper.

For simple3 at SHORT endpoint8, the right cost is1/2 and its
contribution to sigma is1/(2sqrt2). If gamma>0, deleting7 retains
A6, so this is the only left contact and gamma>=1/2. Hence

    delta<=1/2,
    cosh L6>=sqrt12*(1/2+1/(2sqrt2))/sqrt(1/2)
              =sqrt6+sqrt3>11/3.                       (12)

The last inequality follows already from sqrt6>12/5 and
sqrt3>17/10, whose sum41/10 exceeds11/3.

If gamma0, an additional left contact need NOT be absent: these
are the terminal bridges that a B3/F4 sealing argument would miss.
Deleting7 permits at most one contact in A6; while6 is deleted it
lies in A5. The remaining white budget after the right cost1/2 is
less than1/2. The A5 inverse diagonals therefore force exactly a3
at j1 or5 if any left contact occurs: its cost is j(6-j)/12,
and even an endpoint4 would exceed the budget. With no left contact,
or with j1, the selected old5,old6,old7,old8,c path is exact affine
F4. With j5, the full residual and numerator instead give

    delta<=1/12,
    a_6=sqrt12*(5/12+1/(2sqrt2)),
    cosh L6>=5+3sqrt2>9.                              (13)

This estimate retains Q. In particular the A/E finite completions
at the two actual endpoints of that bridge have not been excluded
by an erroneous finite-F4 assertion.

Equations(11)--(13) exhaust the common6/7 inventory. Every surviving
case has an actual edge with cosh L>11/3, and therefore by(3) has
half-tangent square>4/7, a contradiction. B8 pays strictly above4/7.

## Scope

The theorem is geometric and applies to WHOLE B components with
arbitrary old finite spectators at a farthest vertex. Every removal
of a profile used either an exact affine principal, an ACTUAL finite
endpoint, a shadow-forced actual closure, or the explicit Lanner
connectivity argument. No literal reflection group, arithmetic field,
universal finite-index unfolding, or Coxeter property of a metric
facet is imported. It supplies the ordinary rank/14 rate for B5--B8,
not an unrestricted compact or arithmetic dimension bound.
