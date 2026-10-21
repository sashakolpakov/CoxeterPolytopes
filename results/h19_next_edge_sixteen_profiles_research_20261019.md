# The new H19 edge: sixteen finite profiles, twelve after local tests

19 October 2026. Paper analysis and a proposed bounded-search domain
only. No worker was prepared or run for this next edge; no numerical
job, enumeration, Lean process or build ran. Earlier evidence is
unchanged. The28-wall escape is recorded separately in
`h19_twenty_eight_wall_seed_execution_and_escape_20261019.md`.

**Results.** The coordinating agent's16 necessary finite endpoint
profiles are complete. Two short profiles are impossible modulo7,
and two long profiles are impossible modulo the discriminant prime,
leaving11 short and one long necessary profiles. No root in the
remaining12 has been constructed or excluded here. No whole-lattice
nonreflectivity or finite-volume completion follows.

## Fixed actual edge and notation

Let phi²=phi+1, b=4+phi, D=phi*b, F=phi^{-3}=2phi-3, and
omega_j=C_E8^{-1}e_j. The executed walk reached

    V=(1,-phi/2,(phi³/4)*omega0), V²=-F/4,
    d=(0,0,omega1-(7/4)*omega0), d²=7/4, (V,d)=0,
    E=V+d/5, E²=(41-25phi)/50>0.

The fixed eight whites on the edge from V towards E are old0,
old3,...,old8,z. They form mutually orthogonal components

    W=A1(short binary old0)+A6(short chain7--2--3--4--5--6)
      +A1(long z of norm2b).

The same norm-alphabet and finite-angle arguments as the earlier
pinned reports apply. A short/long finite nonright bond is impossible.
Two canonical short roots have contacts0,-1,-phi; two canonical long
roots have contacts0,-b,-bphi, since their product divided by b has
the same allowed twice-cosine alphabet. These statements concern
finite endpoints, where the white-plus-new Gram is positive definite.

## The complete sixteen profiles

A short new root must be orthogonal to the long white z. Write beta
for its positive contact magnitude at old0, and eta for its positive
contacts with the A6. The finite Schur condition is

    beta²/2+eta^T*A6^{-1}*eta < 2.

There is at most one nonzero A6 contact: two unit contacts contribute
at least2, attained by the two endpoints, and the inverse is
nonnegative. At a single position i=1,...,6 along the displayed
chain, the unit contribution is i(7-i)/7, giving6/7,10/7,12/7
in symmetric pairs. A contact phi on A6 is impossible because even
(6/7)phi²>2. Therefore the complete short profiles are:

|Short contacts|Labeled count|
|---|---:|
|No contact|1|
|One A6 contact1, beta=0|6|
|One A6 contact1 at position1,2,5 or6, beta=1|4|
|A6 untouched, beta=1 orphi|2|

The count is13. The middle positions with beta=1 exceed2, and any
A6 contact with beta=phi exceeds2, so no other branch is omitted.

A new canonical long root must be orthogonal to all seven short
whites. Its contact with z is0,-b or-bphi. Dividing both long
normals by sqrt(b) gives the same finite2-node test; each of those
three is positive definite. These are the other three profiles.
Thus the total is13+3=16, with actual white positions retained.

## Short-root equations and two all-height exclusions

For a short candidate q=(x,y,v), put w=C*v. Write alpha=w0,
gamma=w1 and w_j=eta_j at j=2,...,7. Set

    R0=sum_{j>=2}(omega0)_j*eta_j,
    R1=sum_{j>=2}(omega1)_j*eta_j,
    Q0=eta^T*C^{-1}*eta,
    L=4R1-7R0.

Orthogonality to z gives v0=phi³x, hence

    alpha=(phi³x-R0-7gamma)/4,
    y=(beta-phi*x)/2.

Both must be integral, not merely field-valued. Completing the exact
norm equation gives

    B=7gamma+L,
    B²-7F*x²=H,
    H=56-14beta²-28Q0+7R0²+L²
     =28*(2-beta²/2-eta^T*A6^{-1}*eta).             (1)

The block-inverse identity for the free E8 indices0,1 proves the
last equality. It also follows by completing first alpha and then
gamma in the E8 quadratic form. For a single A6 unit contact at
position i, H=56-14beta²-4i(7-i).

If A6 is untouched, R0=R1=Q0=L=0, and the ORIGINAL equation before
completion/division is

    7gamma²-F*x²=8-2beta².                         (2)

The prime7 is inert in O, with residue field F49. Its field norm
sends F=phi^{-3} to-1 in F7. For beta=0 or1, reducing(2) forces
x² to have a nonsquare field norm: respectively -8/F and -6/F have
norm-1 modulo7. A nonzero square cannot have nonsquare norm in F7.
Both profiles are therefore impossible at every height. This is
not inferred from the vacuous completed B²=H reduction alone.

The beta=phi profile is not excluded by that test, and a nonzero
A6 contact has L nonzero modulo7, so gamma occurs linearly in
the original reduced norm equation. No exclusion of those11
remaining short profiles is asserted.

For a root satisfying(1), its edge products are exactly

    (V,q)=-F*x/4,   (d,q)=B/4,
    (E,q)=(B-5F*x)/20.

For x>0 and B>0, H>0 implies B>sqrt(7F)*x>5F*x, since F<7/25.
Thus this branch cuts E and has finite intersection parameter
t=F*x/B, with0<t<1/5 and square-F*H/(4B²). It still requires
every retained wall and angle test; this algebra alone is not a
compatible mirror certificate.

## Long-root equations: only the contact-b branch survives

Write a=0,1,phi for the magnitude of its z contact divided by b.
Orthogonality to the seven short whites and reflection integrality
give every canonical long root in the form

    q=(2u,-phi*u,b*(tau0*omega0+tau1*omega1)),
    u,tau0,tau1 in O,
    -phi*u²+b*(2tau0²+7tau0*tau1+7tau1²)=1,        (3)
    b*(4tau0+7tau1)=2phi³u-a*phi^{-2}.             (4)

Modulo(b)=(D), the residue field is F19 with phi=15 and phi^5=1.
Equation(4) gives u=a/2. Equation(3) then becomes a²=1 modulo19.
For a=0,1,phi the squares are0,1,16. Thus ONLY a=1 survives;
the orthogonal and angle5 long profiles are impossible at all heights.

For a later search in this remaining branch set

    U=2F*u+phi, S=b*tau1.

Completing(3)--(4) gives

    U²-7F*S²=-6bF,
    tau0=(2phi³u-phi^{-2}-7S)/(4b), S/b in O.       (5)

The edge products are(V,q)=-U/4 and(d,q)=7S/4. For u>0, U>0;
the S>0 branch would cut E and have finite parameter U/(7S).
The parameter is less than1/5 by(5), since F<7/25. The other
native integrality, all28 wall and actual endpoint checks remain
necessary. No long candidate is claimed here.

## Proposed bounded domain only: no worker or execution

A later authorized targeted experiment could use the11 short
profiles in(1) with0<x<=600, |x'|<=2, both field-square signs,
then reconstruct alpha,gamma,y and the full root. The single long
branch can use0<u<=300, |u'|<=3/2 in(5), both square-root signs
for S²=(U²+6bF)/(7F), followed by all integrality checks.

These conjugate windows are conservative: at the positive-definite
conjugate embedding the completed native binary norm gives
x'^2<=4/(-D') for short roots and x'^2<=4b'/(-D') for long roots.
Since-D'>2 and b'<4, the short bound2 and the long bound3
(hence |u'|<=3/2) are safe. The physical height600 restriction
would remain only a bounded experiment, not an all-height theorem.

Every accepted root must have all28 native signed products and
Coxeter angles verified, a common strict future timelike center,
a positive cut of E, and a directly checked future timelike
endpoint on the ACTUAL retained edge. The search must retain
both norm classes and must not use an existing-subgroup mirror
orbit. No source or job for this experiment has been prepared.
The present endpoint reduction is16->12, not a chamber completion.
