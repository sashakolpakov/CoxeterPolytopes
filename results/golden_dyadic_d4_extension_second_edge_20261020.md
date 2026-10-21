# The second actual edge of the golden dyadic D4 extension

20 October 2026. Paper proof only; no worker or search was run.
This proves one further actual finite edge, not compactness or reflectivity.

The lattice, complete primitive norm pool, full initial chamber, and first
actual edge are imported from
`results/golden_dyadic_d4_extension_first_edge_20261020.md`, SHA256
`3fc21952a3d119db4370ca46e39aa607666eaea7eaa54d835c3163e7a281035b`.
Thus `phi=(1+sqrt(5))/2`, `O=Z[phi]`, and

    R_n=<-phi>+I_n+D4(O), n>=5,   q in {1,2,4}.

Write the four D4 coordinates last. The first endpoint and its incident
roots are

    Q=(1;0^n;0,0,0,-1),
    B_n + B3(a1,a2,a3) + A1(r),
    a1=e1-e2, a2=e2-e3, a3=2e3,
    r=(2;0^n;0,0,0,-2phi), r^2=4.

All references to finite extensions below concern whole components at an
actual finite vertex of the FULL integral reflection arrangement.

## 1. An actual finite upper and the complete finite inventory

Release a3 and retain `B_n+A2(a1,a2)+A1(r)`. The ray is

    P(t)=(1;0^n;-t,-t,-t,-1), t>0.

The old actual wall `a4=(0;0^n;-1,-1,-1,1)`, of norm4, has
`P(t).a4=3t-1`. It supplies the finite upper `t=1/3`, since

    P(1/3)^2=4/3-phi<0.

The light value is `sqrt(1/(3phi))>1/3`. Hence there is an actual
first new wall no later than a4, and its corner is finite.

For a candidate `w=(a;I;z1,z2,z3,z4)` of normalized norm q put

    c=-w.r=2phi(a+z4),
    b=-w.ai >0 when it contacts one A2 endpoint,
    epsilon=I^2.

The old vertex p is in the chamber closure, so `a>=0`. We allow a0:
the winning wall is an old zero-height wall. A first new wall is not
incident at Q and therefore satisfies `Q.w<0`.

Reflection integrality and the finite native labels give exactly these
individual contacts:

|q|possible positive c at r|possible positive b at A2|
|---|---|---|
|1|none|1, giving label4|
|2|2, giving label4|1 or phi, giving label3 or5|
|4|2 or 2phi, giving label3 or5|2, giving label4|

For q4, `a in2O`, and the four D4 coordinates are congruent modulo2.
These are necessary full-lattice conditions, not merely conditions in a
rational span. In particular c is divisible by2; hypothetical products
1 or phi at the norm4 root r are unavailable to q1.

A new root has at most one contact in A2: two contacts create a cycle.
A connected finite extension of B_n, n>=5, is a B-type diagram. It adds
a norm2 long root with `I=e1`, so epsilon1. It can continue through
one or both nodes of A2, but cannot also reach r, which would introduce
a second4-bond. No q1 or q4 assignment supplies another integral native
extension of this whole B_n. When epsilon0 the complete classification
is the following table; every A2 contact has the two literal endpoint
positions indicated in the last column.

|role|q|epsilon|b|c|multiplicity|
|---|---:|---:|---:|---:|---:|
|isolated|1,2,4|0|0|0|3|
|r only|2|0|0|2|1|
|r only|4|0|0|2 or 2phi|2|
|A2 only, B3|1|0|1|0|2|
|A2 only, A3|2|0|1|0|2|
|A2 only, H3|2|0|phi|0|2|
|A2 only, B3|4|0|2|0|2|
|A2 through new root to r, B4|2|0|1|2|2|
|A2 through new root to r, F4|4|0|2|2|2|
|B_n only|2|1|0|0|1|
|B_n through new root to A2, B_(n+3)|2|1|1|0|2|

There are21 literal profiles. This includes both B4 and F4 bridges.
An H4 bridge would require the unit path `(3,3,5)`: putting label5
at r forces q4, but then the A2 contact cannot have label3, whose
product would be sqrt(2). Putting label5 at A2 forces q2, and the
r contact then has label4; `(3,5,4)` is not finite. The remaining
putative q4 bridge has `(3,4,5)`, also not finite. Thus no H4 type
was silently dropped. These observations, and the classification of
connected finite diagrams containing the retained A2, exhaust the
white-component joins.

## 2. An integral scalar strip with full reconstruction

Put

    j=c/(2phi) in {0,phi-1,1},
    U=a/phi+j in O,      Z=-(z1+z2+z3).

Then

    Q.w=-U,                  a=phi(U-j),
    z4=phi^2*j-phi*U,        U>0, U>=j.                 (1)

For an a1 contact the first three coordinates are `(y-b,y,y)`;
for an a2 contact they are `(y,y,y+b)`. Consequently

    z1^2+z2^2+z3^2=Z^2/3+2b^2/3,
    Z^2-3phi U^2=R,
    R=3q-3epsilon-2b^2-3c^2/4.                         (2)

The intersection parameter is `t=U/Z`. A wall at or before the upper
therefore has

    Z>0, 3U<=Z,
    U^2<=R/[3(3-phi)],       U'^2<=phi R'/3.           (3)

The second inequality uses the other real embedding in(2); all the R
in the table below are totally positive. No old-wall acuteness was
used to make this necessary domain smaller.

For clarity the exact omitted glue, to be tested after this strip, is

    b0: Z in3O;
    a1: Z=b mod3, y=(b-Z)/3;
    a2: Z=-b mod3, y=(-Z-b)/3;
    D4: z4-Z in2O.                                      (4)

Together with the coordinates in(1) and the q4 conditions above, these
recover the full actual lattice vector. There is no replacement of
a whole component by an embedded subcore.

The R-values of the complete inventory are

|role|R|
|---|---|
|isolated q1,q2,q4|3,6,12|
|r only q2,c2|3|
|r only q4,c2|9|
|r only q4,c2phi|9-3phi|
|A2 only q1,b1|1|
|A2 only q2,b1|4|
|A2 only q2,bphi|4-2phi|
|A2 only q4,b2|4|
|either A2/r bridge|1|
|B_n only|3|
|B_n/A2 bridge|1|

## 3. Complete tiny coefficient census

Write `U=A+Bphi`, so `B=(U-U')/sqrt(5)` is an integer. The exact
bounds(3) give

|R|complete positive U list|
|---|---|
|1,3,4,4-2phi|empty|
|6|1,phi-1|
|9|1,phi-1|
|12|1,phi-1,phi|
|9-3phi|1,phi-1|

Here is a short complete coefficient check, rather than a height
cutoff. The elementary inequalities
`8/5<phi<81/50` and `11/5<sqrt(5)<56/25` suffice throughout.

* R1 has `U<1/2`, `|U'|<3/4`, so B0 and no positive integer U.
* R3 has `U<6/7`, `|U'|<13/10`. Their sum is less than sqrt(5),
  again forcing B0 and no positive integer.
* R4 has `U<1`, `|U'|<3/2`, hence B0 or1. The only remaining
  candidate phi-1 has conjugate absolute value phi>3/2.
* R=4-2phi has `U<1/2`, `|U'|<17/10`; their sum is less than
  sqrt(5), again leaving no positive integer.
* R6 has `U<5/4`, `|U'|<9/5`, hence B0 or1; these give only
  1 and phi-1. Both satisfy the exact strip.
* R9 has `U<3/2`, `|U'|<sqrt(5)`, hence B0 or1, with the same
  two values.
* R12 has `U<7/4`, `|U'|<13/5`, hence B=-1,0,1. For B=-1,
  the possible values 2-phi and3-phi have conjugate absolute value
  at least phi^2, whereas the exact conjugate bound is `2sqrt(phi)`.
  Since `phi^4-4phi=2-phi>0`, both fail. B0 or1 gives the stated
  three values, which satisfy the exact strip.
* R=9-3phi has `U<=1`, `|U'|<=sqrt(3phi+1)<5/2`, so again
  B=-1,0,1. The only possible B=-1 value is2-phi; it fails since
  `phi^4=3phi+2>3phi+1`. B0 or1 gives the stated two values.

For the nonempty R6 strip, U1 gives

    Z^2=6+3phi,             N(Z^2)=45,

which is not a square in K. At U=phi-1 it gives `Z^2=3phi^2`;
this is not a square because3 is inert in K. For R9 the two values
give `Z^2=9+3phi` and `6+3phi`, of rational norms99 and45.
These are nonsquares as well. Thus both R6 and R9 disappear before
any glue relaxation could introduce a candidate.

R12 occurs only for isolated q4. Here j0, and(1) gives respectively
`a=phi,1,phi^2`, all units. None lies in2O, so every R12 possibility
violates the full q4 reflection condition.

Finally R=9-3phi occurs only for q4 with c=2phi, hence j1. The
old p-side condition `a>=0`, or `U>=1`, leaves only U1. Equations
(1)--(2) then give

    a=0, z4=1, Z=3, (z1,z2,z3)=(-1,-1,-1),

exactly a4. It satisfies all full-lattice conditions and every old
chamber inequality at the displayed intersection. Therefore a4 is
the unique first new simple wall; no simultaneous ambiguity remains.

## 4. Actual conclusion and the remaining gap

The next actual vertex is

    S=P(1/3),      whole diagram B_n+A2+I2(5),

where the I2(5) pair is r,a4. This is a complete second edge of the
full arrangement, uniformly for n>=5, and in particular for the
dimension9 candidate R_5.

It proves neither compactness nor an infinite chamber period. In
particular the reflection interchanging r and a4 would have normal
`r-a4`, whose norm is `4(phi+2)` and whose vector is primitive in
R_n. Its norm has a prime above5, outside the complete primitive
norm pool. Thus the formal I2(5) reversal is not by itself an
integral arrangement symmetry. Further outgoing-edge or module
geometry is required; no unbounded prefix or new example is claimed.
