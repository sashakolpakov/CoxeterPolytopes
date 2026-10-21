# The triple-H3 corner has a full compact tetrahedral face

20 October 2026. Exact paper proof only; no new worker, search or
compact dimension-nine example. The sixteen-wall prefix is already
closed negatively by its spacelike witness. This note identifies the
actual face suggested by that witness, rather than extending a prefix.

Import the frozen four-edge and full M6 chamber packet
`results/golden_dyadic_d4_research_20261020_freeze.json` and use the
coordinates of
`results/golden_dyadic_sixteen_wall_caps_20261020.md`,
SHA256 `c288b5de619471772c8381d98afa3996ac8851b4d2d146ad520d6474ee8ba92b`.
In particular the actual H3 of primitive norm4 roots(a3,s,r) is fully
sealed, and its entire actual face is the full native M6 chamber.

## 1. A different actual corner and its exact perpendicular lattice

In M6 take

    M=(1;phi^-1,phi^-1,phi^-1,phi^-2,phi^-2,phi^-2),
    M^2=21-13phi=-phi^-7.

The active native M6 roots split into two orthogonal components

    E=(u,beta1,beta2),     F=(w,beta5,beta4),

both H3 of primitive norm2, with bonds5,3. All other M6 wall values
are strict: beta3=-phi^-3 and beta6=h=-phi^-2. Since the full M6
face is already identified, this is an ACTUAL ambient finite vertex
with complete diagram `H3_2+H3_2+H3_4`, not a Gram-only proposal.

Put

    T0=phi^4*M=(phi^4;(phi^3)^3,(phi^2)^3), T0^2=-phi.

This is an integral unit-norm time vector, so

    M6=O*T0 orthogonal-sum N6

integrally, where N6 is positive definite and unimodular. The E and F
root spans lie in N6. Their determinant is `2phi^-2`. Each span is
saturated in its K-space: a nontrivial integral overlattice would have
index ideal whose square divides(2), which is impossible. Hence
`N6 intersection K E=E`, and similarly for F.

Here is the full seal, including nonparabolic finite overgroups. A
mirror crossing `(E+F)^perp` generates with E,F a finite reflection
group at its timelike intersection. The finite classification allows
only `H4+H3` for a new nonorthogonal component: a component containing
a5-bond is H3 or H4, and a rank4 H4 cannot contain the two independent
rank3 components together. In the H4 case every H4 root is orthogonal
to F. Its primitive length is2 by integral reflection conjugacy with
the original E roots. The H4 simple-root lattice is unimodular.

But any ambient integral vector orthogonal to F has its E projection
in E. Indeed remove its integral T0 and old norm4-H3 components; its
remaining M6 component lies in `N6 intersection K E=E`. Thus projection
splits E integrally from the proposed unimodular H4 root lattice.
Its nonunit determinant `2phi^-2` cannot divide a unit determinant.
This excludes H4, even if E is not assumed parabolic in that H4.
Every crossing mirror is therefore orthogonal to BOTH E,F.

The exact perpendicular lattice is consequently

    T3=O*T0 orthogonal-sum P3.

It splits integrally from R5, and the full intrinsic reflector
arrangement equals the induced full arrangement. This is an actual
three-dimensional face. This seal is special to the whole E+F context;
it does not assert that either norm2 H3 is sealed in arbitrary contexts.

Write `wP=a2-V0=-2omega1`. An exact O-basis of P3 is(a3,s,wP), with

    Gram=[[4,-2,-2],[-2,4,0],[-2,0,phi+2]],
    determinant=4phi^4.

The determinant ideal equals that of P3, so no omitted glue remains.
In particular the vectors of P3 orthogonal to s,r are exactly O*wP.
The full primitive norm pool stays1,2,4 by the integral orthogonal
split and the imported norm pool of R5.

## 2. One complete actual edge, with no height search

Release a3, retaining the two norm2 H3 components and the norm4
`I2(5)(s,r)`. The explicit integral norm2 upper is

    b=T0+wP=T0-2omega1.

The ray in T3 is `T0-tau*omega1`. Since
`omega1^2=(phi+2)/4`, the upper meets it at

    tau_b=2phi/(phi+2)=2/sqrt(5),

strictly before light. This is an upper, not yet an actual-wall claim.
Its products with every wall of the finite sixteen-wall list are
native or dotted: beta3=-phi, beta6=h=-phi^2, a3=-2,
a2=A(a2)=-4phi-1, W=-phi, Z=-1; all other products vanish.
It cuts the old escape E strictly, with `b.E=phi^-2>0`.

Let n be ANY first new simple wall on this actual edge. By the full
seal it has the form

    n=a*T0+p, a in O, a>0, p in P3, n^2=q in{1,2,4}.

Let L=p.wP. Contact before or at b gives L>0 and
`L/a>=phi+2`. Positive Cauchy--Schwarz in P3 gives

    p^2>=L^2/(phi+2),  hence a^2<=q/2.               (1)

At the other embedding the entire form is positive definite, so

    a'^2<=q*phi.                                    (2)

These inequalities are a COMPLETE finite interval bound, not a
selected root-height cutoff. For a=A+Bphi,

    Tr(a^2)=2A^2+2AB+3B^2.

If q1, (1),(2) give trace at most `1/2+phi<3`. This permits only
a=0,+1,-1; none satisfies `0<a^2<=1/2`. If q4, reflection integrality
gives `n.T0 in2O`, hence a in2O. Applying the same q1 bounds to a/2
excludes it. Thus q2 is the only remaining class.

For q2, the trace is at most `1+2phi<5`. The nonzero possible squares
come from `a=+/-1,+/-phi,+/-phi^-1`; the physical bound `0<a<=1`
leaves only

    a=1 or phi^-1.

A norm2 new root can meet a norm4 old root at a finite native angle
only by a4-bond, with product-2. Either attachment to the retained
I2(5), or both attachments, is nonfinite (the3-node4,5 chain already
has negative determinant). Thus n is orthogonal to BOTH s,r. The
exact P3 basis then gives `p=m*wP`, m in O. Its norm equation is

    (phi+2)m^2=q+phi*a^2.

At `a=phi^-1` this says `m^2=phi^2/(phi+2)`, which is not integral:
its field norm is1/5. At a1 it gives m^2=1; positive L forces m1.
Consequently b is the UNIQUE first simple wall. No unexamined
simultaneous positive root or bounded-height remainder is omitted.

## 3. The full lower face is compact, so it is not an obstruction

The four ACTUAL intrinsic walls(b,a3,s,r) form the chain

    [4,3,5].

Their full Gram determinant is `-16/phi`. Each3-wall subdiagram is
finite: respectively H3, `I2(5)+A1`, `I2(4)+A1`, and B3. Thus they
bound a compact hyperbolic Coxeter tetrahedron. Equivalently the
four positive definite3-wall completions give its four finite
vertices. All are feasible for the full chamber by the positive
acute-projection argument for actual walls, and their convex hull
is the whole tetrahedron. Therefore it is the FULL T3 chamber.

This proves a small exact compact face of the open dimension-nine
host. It prevents using this face for a noncompactness obstruction
or an infinite-order full-chamber symmetry. It neither proves nor
disproves compactness of R5. A genuinely new global step must leave
this fully identified compact face rather than add an unbounded
sequence of prefix caps.
