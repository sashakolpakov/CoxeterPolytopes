# Independent written review: the L59 D7+c first edge

20 October 2026. **PASS; no mathematical flaw found.** I independently
derived the finite extensions and their norm equations before reading
the complete interval note, then read its entire text and checked both
ten-entry norm tables by hand. No scientific execution, imported worker,
census, Lean or build was run. Actual checker execution is not claimed.

Reviewed interval source:

* `compact_h9_golden59_d7_edge_interval_20261020.md`
* SHA256 `e4715a78f22628bf7a27051ae6588adae7803f191a2a08f476d0a02fa2223f52`

Inherited route interface: SHA256
`725d8775cebec51d8c18b4b0a0d3be0255e8b1472621adb4d3cefa656adb429e`.
The separate period proof is SHA256
`699e9922fc0cc3ca6f16d5619e02ab345fbcb39ed25ddb355d32a4120ef14f9d`.
I authored that period proof, so this is an independent census review,
not an additional independent review of my own isometry derivation.

## Finite inventory and reconstruction

The retained diagram is D7 with arms1,1,4, plus the short isolated c.
An attached short root either makes D8 at the long end, E8 at one of
the two fork ends, or D9 by bridging the long end to c. Other attachments
produce a cycle, an extra fork with an affine D-subdiagram, or a
nonfinite arm triple. A5-bond cannot occur in a connected rank>=7
finite diagram. If D7 is untouched, the c contact is0,1 orphi. The long
class is orthogonal to all whites by the corrected whole-arrangement
isolation argument. This gives exactly the eight literal profiles;
no suppressed spectator contact or extra finite rank-nine component
was found.

I rederived m7=2,m1=5,m2=7 and white projection costs1,7/4,7/4. Writing
X=2 times the h coefficient gives precisely

    Z^2-2delta X^2 = 8,6,6-2phi,8b,4,1,1,2

in the stated order, and tau=delta X/Z. Modulo b, the four excluded
residues are8,6,56,2. Quadratic reciprocity gives2 nonsquare,3 square
and-1 nonsquare in F59, so all four exclusions are all-height. They
are not finite-box observations.

For survivors k=0, h_y=phi^3 is a unit and integral binary coordinates
force X in2O. The E8 spatial congruences follow from its unimodular
Cartan matrix, not an arbitrary choice of denominators. The long row's
spatial coefficient is b times an integer, and its divided equation
indeed gives X^2=45 modulo59, hence X=+/-24 modulo b.

## Complete interval and both finite lists

The proposed long upper is integral and reflective: Gh=(-delta,delta),
its E8 covector is b(9+14phi)e0, and its norm is2b. The useful identity
delta*t^2-2z^2=-b holds for t=38+61phi,z=133+213phi. Its contact is
strictly timelike. Monotonicity of the positive Pell branch gives exactly
8b X^2<=R Xv^2, with the direction of the inequality as printed.
The conjugate bounds and the physical bounds25,50,274 are conservative
in the safe direction.

For D8 the substitution X=2Y,Z=2W gives
W^2=1+2delta Y^2. For E8 the same square test uses Y=X directly, so the
common larger list is valid even before parity/congruence filtering.
The exact inequalities imply0<=d<=12 for Y=a+dphi. To check its ten
members without trusting a census, use the rational bracket
1.618<phi<1.619. The allowed conjugate radius is between.45 and.46.
Among d=0,...,12, nearest integers to(phi-1)d retain exactly
d=1,2,3,5,6,7,8,9,10,11; d=4 misses that radius, d=12 violates Y<25,
and d=0 gives no positive entry. All retained remainders have absolute
value<.45, and the excluded d4 remainder has absolute value>.46.
These are exactly the ten printed Y values.

For the long row the exact equation gives
A^2=(-4+14phi)t^2+1-phi, with0<t<=9+14phi and
2b' t'^2<=1. Its conjugate radius lies between.35 and.36. For
0<=d<=16 the retained d values are2,3,5,6,7,8,10,11,13,14; the
remaining positive candidates either miss the radius or exceed the
physical upper. The d14 remainder is less than.35, so it is genuinely
included. This independently reproduces the second ten-entry list.

For Y=a+dphi, put n=a^2+ad-d^2. I checked the first norm table through

    N(1+2delta Y^2)
      =1+44(a^2+d^2)+82(2ad+d^2)-236n^2.

I checked the second table through

    N((-4+14phi)t^2+1-phi)
      =-1+38(a^2+d^2)+44(2ad+d^2)-236n^2.

Every printed integer agrees. All are nonsquares except1681=41^2
at t=9+14phi; the corresponding radicand is exactly
(38+61phi)^2. Both possible signs of the field square root are accounted
for, and X>0 chooses the positive one. Thus the full surviving interval
has only the primitive long upper, independently of cap compatibility.

## Actual first wall and global-proof boundary

The existing actual edge starts at the full E8+c vertex. Its first
intersection with the locally finite FULL integral arrangement occurs
no later than the proved finite upper. At that point the chamber has
a finite simple vertex, so its new simple wall must belong to the
exhausted spherical inventory. Therefore the unique surviving long
root is the actual first simple wall, with full germ D7+c+v. There is
no simultaneous-candidate ambiguity. This argument needs no prior
compactness or finite-volume claim about the whole chamber.

The postselection requirements (native products with the sixteen actual
old walls, full reflection integrality and positive future apex) remain
explicit checker interfaces, not filters used to manufacture uniqueness.
The displayed H has the required norm and orthogonality by the exact
delta*t^2-2z^2=-b identity. The written first-wall argument supplies
the missing genuine-germ premise of the separate infinite-period proof.

This review supports the written literal-host exclusion when combined
with the independently reviewed period proof. It does not record an
executed exact certificate, a new compact9 example, a new dimension
bound, or a statement about the full commensurability class.
