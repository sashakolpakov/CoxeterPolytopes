# Independent review of the H19 full edge and period

20 October2026. **PASS: no mathematical flaw found.** I read the complete
paper and independently rederived its interval, chamber maps and trace.
No scientific program, import, search, checker, Lean or build was run.

Reviewed source `compact_h9_discriminant19_d7_edge_period_20261020.md`,
SHA256 `ecc61a2dcff6df1678cc94533aa201c70cd0b262e26deab5ef85f5d64ebaf824`.
I also read the full inherited mixed-norm report70fa0107, polar-cap
report792634e9, the corrected whole-arrangement prime59 isolation
proof725d8775, the earlier facet-reversal argument1d9a7124 and the full
`bugaenko_vinberg_transcript.md`113e5469. Their exact archive checks are
inherited, not independently executed in this review.

## Full arrangement and the first transverse edge

The prime19 isolation calculation uses the complete finite-angle list,
including the two nonfundamental golden values2-phi and3-phi. The listed
norms of those values and of their differences from4 are correct and
exclude19 divisibility. For two primitive long roots their nonzero
reductions occupy the same one-dimensional radical modulo b; the norm
and pairing calculations modulo b squared give a squared-angle residue4.
Thus distinct long mirrors cannot meet, and a short mirror meeting one
is perpendicular. Anisotropy also excludes parallel distinct mirrors.

This upgrades the r-facet to the complete restricted norm-two arrangement
of M=<-2phi>+E8. I checked that the inherited Vinberg transcript supplies
the needed full-chamber input: its separate M norm-pool argument excludes
primitive norm4, its ordered enumeration includes every height through
the last wall, and its independent41-vertex compactness certificate
allows stopping. The proof does not replace that input with the false
assertion that an arbitrary compact root prefix is complete.

The binary roots r,c are an integral basis. A root crossing the open
r-to-c segment has coefficients lambda,mu of opposite signs. In the
short case the physical norm gives lambda squared<=1/b; completing the
conjugate square gives lambda prime squared<=4/[b'(4-b')]. Their
product is below1, contradicting the integral norm of nonzero lambda.
For a long root, pairing divisibility gives mu=b gamma; division of the
norm by b gives exactly the same obstruction for gamma. I checked both
completions of squares, the factor19 and4-b'=phi-1. Consequently c is
first in the full arrangement and p has the stated actual E8+A1 germ.

## Complete D7+c interval

Finite extensions of the D7 with arms1,1,4 give the eight placements in
the source. The only short connections to c are0,1,phi. The connected
D7 extensions are D8, the two E8 forks and the D9 long-arm bridge;
long roots are orthogonal to every white. No contact with a spectator
or exceptional finite rank-nine component was silently discarded.

For the native binary coordinates x=X,2y=k-phi X and the displayed
weight projection I independently obtain

    Z squared-2D X squared=4q-2k squared-4q_white=R.

The E8 residual norms are1 at node6 and7/4 at nodes1,7. Thus the eight
R values, including the long value8b, are correct. Modulo b, phi=15;
the residues8,14,2 are nonsquares in F19. These eliminate precisely
the isolated short, golden c-bond and D9 bridge, without a height cut.

The polar z is a primitive integral long reflector: its y coordinate is
a unit and its covector is divisible by b. I obtain its norm2b,
tau=phi cubed/2 and endpoint square3-2phi<0. Comparing positive collision
ratios gives K=phi to the eighth/(2b)=(44+71phi)/38. Conjugate positivity
gives the separate X-prime bound. All rational coarse bounds in the
source follow from8/5<phi<13/8.

For k=0, X=2Y follows from exact binary integrality. The E8 strips force
Y to be an integer with conjugate absolute value below1/4, impossible.
For D8 they allow only the noninteger candidate Y=phi squared, which
violates the physical bound. For A2 the complete coefficient strip and
the congruence phi X=1 modulo2 force X=phi squared; the resulting square
would have rational norm704, strictly between26 squared and27 squared.
No spatial-divisibility filter is needed for these exclusions.

For the long row, the primitive covector condition and E8 unimodularity
give t=Z/(4b) integral. With X=2A its exact equation is

    2b t squared-phi A squared=1.

The interval gives0<t<=phi squared and the other embedding gives
|t prime|<2/5. The full coefficient strip then forces t=phi squared,
and positive X selects A=phi to the fourth. This reconstructs z only.
Every finite profile is therefore covered, and local finiteness plus
the finite upper proves a unique first simple wall. The later31-wall
prefix and its escape are not used. The new full germ is D7+c+z.

## The two complete chamber maps

The displayed P has norm-2phi and covector(-phi to the fifth,1;
phi squared e0). Its point inversion is integral because phi is a unit.
I directly obtain J_P r=z and J_P c=g. For the actual comparison vertex
B, the walls D7+r+g are positive definite and all are actual walls:
g belongs to the complete compact r-facet, not merely to a seed list.
The acute positive-definite projection argument makes their intersection
feasible and strict on all other walls. This supplies the full ambient
finite germ, also when C is not yet known to have finite volume.

I checked the stated longest-D7 action, including its e0 column and its
fixing omega0. Thus H=J_P w0 is integral, future-preserving and swaps the
two outward D7 forks. It maps z to r and c to g. The identity
HV=P-V=(phi-1)B has positive scalar, so H maps complete inward actual
germs of the same chamber and hence preserves C.

The A columns preserve the form and fix r. The extension column on e_y
is exactly(r-Ae_x)/2 and is integral as displayed; no index-two extension
is assumed. I independently checked A squared=I and the oriented facet
permutation, including -e7 to k and g to old9. Since its y coordinate is
unchanged, A preserves the occupied side of the complete r-facet.
Local full-arrangement facet adjacency therefore gives A(C)=C.

## Trace and literal-host conclusion

For A w0 I obtain exactly the five nonzero diagonal contributions

    9+16phi, 1, -17-30phi, 12+20phi, -1,

whose sum is4+6phi. With theta=(2,4,6,5,4,3,2,3), I independently find

    AP=(8+13phi,-phi;phi squared theta+phi to the sixth omega7).

Its first spatial coordinate is27+42phi and its pairing with P is
-20-34phi. The rank-one point-inversion trace formula consequently gives
trace(AH)=10+14phi>10. This proves infinite order; it is not inferred
merely from two distinct involutions.

Both factors preserve the same full chamber. The standard free action
of the integral reflection group on its open chambers makes the powers
of AH distinct cosets in the admissible arithmetic lattice O^+(L).
Thus the full reflection group has infinite index and infinite covolume.
A compact Coxeter reflection group preserving this literal L would be
a lattice subgroup of that discrete overgroup, which is impossible.

The conclusion is nonreflectivity of this literal discriminant19 H9
lattice. It is not a uniform dimension improvement, a compact example,
or an exclusion of its entire similarity or commensurability class.
The standard finite spherical classification, admissible arrangement
local finiteness, actual-germ continuation and arithmetic index/covolume
facts remain explicit written inputs. Exact archive and future worker
calculations supplement these layers and do not replace them.
