# A common center/long-neighbor opposite in a D component

20 October 2026. Complete paper proof, no scientific computation,
optimizer, build or formalization. Other finite components are arbitrary.

## Statement and notation

At an actual vertex of a compact Coxeter polytope suppose its WHOLE
finite diagram has a D_r component, r>=7. Use standard unit roots
`(e_i-e_(i+1))/sqrt(2)` and `(e_(r-1)+e_r)/sqrt(2)`. Write

    o=r-2,     a=r-3,     b=r-1,     e=r,

where o is the fork center, a its long-arm neighbor, and b,e the two
fork leaves. If ONE actual outside wall c is the actual opposite of
BOTH directions o,a, its COMPLETE old contacts, up to exchanging b,e,
are exactly one of the following possibilities:

    A: contact4 at a and3 at o;
    B: an isolated contact5 at b;
    C: contact3 at a and4 at b.                         (1)

Here “complete” includes every old orthogonal spectator component:
all unlisted products are zero. The conclusion is a necessary
classification; it does not assert all listed profiles occur in every
rank or in a compact polytope.

Each possibility forces the D_r component covering energy to exceed
3/4. Consequently a D7 component has covering energy greater than1/2
at any actual compact Coxeter vertex, with arbitrary spectators. At a
farthest vertex this is a component estimate in the additive energy
less than1, not a new ambient dimension ceiling.

The covering-energy and exact actual-opposite conventions are those of
`compact_universal_e6_e7_e8_component_energy_20261020.md`, SHA256
`acf261d97f6eb1a11a64502e17e079ef7f418b9746fbf85dbdf6915750fe9a18`.
Its sharper pair inequality and explicitly external classical input
that Lanner diagrams have order at most5 are retained. Actual affine
principal diagrams are impossible in a compact polytope. Finite
connected Coxeter diagrams and their A/B/D/E/H continuations are used
in the elementary small-chain arguments below.

## Complete contact inventory from the two actual endpoints

Deleting a leaves the right-hand A3 path b-o-e, the left A_(r-4),
and all spectators. Thus the possible contacts of c to b,o,e are:
none; a single3 at o; or one contact3,4,5 at a fork. There cannot be
two: a middle/endpoint pair gives a retained nonpositive triangle,
and two endpoints give a retained nonpositive four-cycle. A4 or larger
contact to o gives a retained nonpositive star; a6 or larger at a fork
contains an affine6--3 path. This exhausts the right-hand contacts.

Deleting o leaves the whole left A_a chain1,...,a, the isolated old
forks, and the spectators. Positive definiteness permits at most ONE
contact from c to that left chain: two would give a nonpositive cycle.
All subsequent exclusions use actual finite replacements missing a
or o, or exact affine diagrams. A nonpositive diagram in the full
star is not treated as forbidden without that additional reason.

**Right contact3 at o.** If c.a=0, then c,o,a,b,e give an exact
affine D4 star. Hence the contact at a is positive. In the replacement
missing o it is at the endpoint of an A_a chain with a>=4, so it is
3 or4. Contact3 gives the exact affine triangle c,a,o; hence it is4.
Any second left contact is already excluded by the retained cycle.
No other white contact is possible: in the replacement missing o,
extending the 4-end c produces the nonpositive five-node path
`3,4,3,3` (affine F4 at the smallest native new bond). The required
old three-node chain exists since a>=4. Thus (A) is complete.

**Right contact5 at b.** In the replacement missing a, c,b,o,e form
H4. No finite connected Coxeter diagram extends H4, so c has no other
WHITE contact there. A positive contact at the released a would, in
the replacement missing o, give a chain at least as strong as
`3,3,3,5`, which is nonpositive. Its raw norm2 determinant is
`6-4phi<0`. Thus this contact is zero as well, giving (B).

**Right contact4 at b.** In the replacement missing a, c,b,o,e form
B4 with c at its 4-end. Any additional white contact to c gives the
nonpositive path `3,4,3,3`; thus there are none. If c.a=0, the old
center with arms a,e,b-c gives an exact affine terminal-4 fork. A
positive contact at a must be3: a4 or larger gives a retained
nonpositive4--4 path with b after deleting o. This gives (C).

**Right contact3 at b.** A positive contact at a is impossible:
contact3 gives the exact affine four-cycle c,a,o,b; contact4 or
larger gives, after deleting o, a path at least as strong as affine
F4 `3,3,4,3`. Any contact at another left node j is also impossible.
A3 contact makes the exact all3 cycle through c,b,o and the old path
to j. An internal4 or larger has a nonpositive retained star; an
end4 at j=1 gives the retained affine F4 path old3,old2,old1,c,b.
A5 or larger end contact is already nonpositive in the retained
long chain. Consequently c has no left contact.

If c contacts a spectator node x, the finite replacement missing a
contains the chain x,c,b,o,e. Its new bond is3 or4: a5 or larger
makes the retained path5,3,3,3 nonpositive. Contact3 makes the full
old diagram contain an exact affine E7 star with arms1,3,3, using
e, the arm b-c-x, and three nodes of the old long arm. Contact4
makes the exact affine terminal-4 fork with arms1,1,3, using e,a
and b-c-x. Thus there is no spectator contact either. With no other
contact c would be isolated3 at b: for r=7 its own projection cost
is r/8<1, so its height square would be negative; for r>=8 the full
star contains the exact affine E8 tree with arms1,2,5. This excludes
the entire right-contact3 case.

**No right contact.** There is no contact to o,b,e. If there is no
left contact either, the directional derivatives at o,a vanish, so
c cannot close them. Otherwise let j be the unique left contact.
The finite replacement missing o permits a3 contact, or an end4.
A3 at an internal old node j>=2 creates the exact affine two-fork
D diagram: one fork has leaves c,j-1, and the other has leaves b,e,
with the old path from j to o between them. An end4 at j=a gives
the exact affine4,3,3 star at a; an end4 at j=1 gives the exact
affine terminal-4 D-extension along the long arm. All are forbidden.

The sole remaining case is a3 at j=1, possibly with arbitrary
spectator contacts. Choose a minimal nonpositive principal subdiagram
of the full old vertex plus c. It is Lanner, since the polytope is
compact. It contains c and o: otherwise it survives deletion of o,
contradicting that actual finite replacement. It is connected. But
the unique old path from c to o has r-2>=5 edges. Spectator nodes
attach only through c and cannot shorten it. Thus a connected
subdiagram containing c,o has at least6 nodes, contradicting the
Lanner order bound. This completes every left/right/spectator case
and proves (1).

Whenever a larger bond or extra contact was said to worsen a retained
obstruction, the base connected nonpositive Gram has a positive
nonpositive quadratic witness. Decreasing its off-diagonal entries
cannot restore positive definiteness. The exact affine exclusions
above instead have all their displayed internal products fixed.

## Actual edge pressure for all three profiles

Put k=r-2=a+1 and use the unit inverse-D entries

    C_aa=2a, C_oo=2k, C_ao=2a,
    C_bb=r/2, C_ab=a, C_ob=k.

For (A), the positive white residual at the actual o endpoint is
delta_A=1/k and the directional derivative is `(sqrt(2)*a+k)/sqrt(2k)`.
For (C), the residual is delta_C=1/(2k) and the derivative is
`(a+k/sqrt(2))/sqrt(2k)`. Thus BOTH give exactly

    cosh L_o=a+k/sqrt(2)>7,     r>=7.

Their actual half-tangent square is `(cosh L_o-1)/(cosh L_o+1)>3/4`.

For (B), put beta=cos(pi/5), so beta^2=(3+sqrt(5))/8. Its height
square is r*beta^2/2-1 and its a-direction derivative square is
beta^2*(r-3)/2. Consequently

    tanh^2 L_a=1-(9-4sqrt(5))/(r-3).

This is increasing in r. At r=7 its half-tangent square is EXACTLY
`4/sqrt(5)-1>3/4`, with strict seal `256>245`. Hence every profile
in (1) forces component covering energy greater than3/4. These
formulas use full and half tangents separately.

## D7 energy greater than one half

The center/long-neighbor correlation in D7 is k0=sqrt(4/5). If its
component covering energy were at most e=1/2, the sharper pair bound
from the pinned source would give

    tanh L_o*tanh L_a <=16(1+k0)/(5+k0)^2 < k0.

Indeed the cleared difference is `k^3+10k^2+9k-16`, increasing for
k>=0; at k0 it is `(49/5)k0-8>0`, sealed by `9604>8000`. Distinct
actual opposites require the reverse inequality, so o,a must have
the SAME actual opposite. The complete classification and pressure
above would then force energy greater than3/4, a contradiction.

Thus E_D7>1/2 with arbitrary spectators. For larger D components,
distinct o,a opposites remain an unresolved branch. No ambient
dimension ceiling or all-D energy rate is asserted.
