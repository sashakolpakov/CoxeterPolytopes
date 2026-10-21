# Universal order11 B/C zero transfers

The order11 family costs can become A99/2, B53, C52 and D53. In affine
B10, delete the two remaining terminal/fork halves. In affine C10,
delete the two distance9 pairs next to its already deleted terminal
pair. Add one explicit finite upgrade: a whole A1 against the doubled
terminal of a separate whole B10 now pays1 instead of1/2. Its actual
added cost is1/2 per such component pair.

Every affected quadrilateral still receives at least1 and triangles
retain at least2. The new rank10 ideal envelope is53/11. This is a local
payment theorem; a simultaneous finite count and global scalar comparison
remain separate. No new dimension bound is asserted.

The [source](../scripts/fresh_order11_bc_zero_transfer_20261018.py) and
[certificate](fresh_order11_bc_zero_transfer_20261018.json) build on the
[coupled B/C/D half theorem](fresh_order11_bcd_half_pair_transfer_20261018.md),
the independent [C endpoint-zero theorem](fresh_order11_c_endpoint_zero_transfer_20261018.md),
and the [completed dimension61 rule](fresh_rank15_rank16_dimension61_20261018.md).
They also bind the new [order16 D105 theorem](fresh_order16_d_distance9_half_transfer_20261018.md).
All root norms are2, with native simple magnitude1 and doubled magnitude
sqrt2. No arbitrary rank or label cutoff is used.

## Actual-star identities and direct alternatives

At an affected ideal corner, let a,b be its marks and c,d their opposite
black roots. The actual varying affine factor F has positive kernel nu;
the remaining supporting whites form a positive-definite block Z,
orthogonal to F. Every actual edge is strictly positive definite and
every actual corner positive semidefinite. These are the inherited
actual factor-fiber facts. Write x,y for the nonnegative contacts of
c,d to Z and h,k for their positive heights against nu.

An adjacent black contact gives a distinct marked path1; a common
spectator contact gives a distinct marked path2. Both remain full units
under every rule in this packet. Remove those branches. Adjacent black
products now vanish and x,y have disjoint supports.

After eliminating F minus its marks, the scaled old-mark Gram is
R^(-1)[[1,-1],[-1,1]]. For each opposite root let l be its moment and
c0 its residual constant. Put g=2l/h+c0/h². Eliminating Z from the
actual three-dimensional face quotient gives the determinant identity

    0 = 1−(g_c+g_d)/R
          +(x/h−y/k)^T G_Z^(-1)(x/h−y/k)/R.             (1)

The final form is nonnegative. The source checks this symbolic identity
at R9 and R10, as well as the physical profile and cross identities.
It rejects floating-point constants in every saved symbolic value.

## C: the near-terminal distance9 orbit

Number the affine C10 path0,...,10, with doubled bonds0–1 and9–10.
The already deleted pair is(0,10). The new zero orbit is(0,9),(1,10),
invariant under reflection. Consider a=0,b=9; the reflected case is
identical. Kernel coefficients are1 on internal nodes and1/sqrt2 at
the two terminal nodes. The exact resistance in(1) isR=9.

Deleting a leaves the finite B10 path1,...,10. Its complete PSD
one-root extensions have no contact, a simple or doubled contact at
its near endpoint1, or a simple contact at node2. The four heavy-end
cost bounds are respectively4,3,2,4, including every possible contact
to the omitted opposite mark. These are the same complete four
profiles as in the endpoint-zero theorem; the source rechecks their
physical identities at the new resistance9.

Deleting b leaves finite B9 on0,...,8 and an isolated white10. For the
B9 block the norm-two inverse diagonal is1 at8,2 at7, and greater
than2 everywhere else. An extension having two contacts inside B9
would contain a cycle. A connected finite diagram is a tree; the only
affine cycle type is A, which cannot contain the old doubled bond.
Exceptional finite/affine types are too small to contain this ten-node
extension. Thus at most one B9 contact is possible. The Schur condition
allows only simple/doubled at8 or simple at7, including equality.

Let q be the magnitude of the contact with old white10, and t>=0 the
opposite product with the omitted mark b. Strict positivity of the
actual edge containing10,d gives q<2; a nonzero native q is at least1.
The PSD corner containing the complete B9+A1 block further bounds the
sum of its projection costs. The complete necessary profiles are:

|B9 contact|Additional q|Height k|Moment|c0|Maximum g|
|---|---|---|---:|---|---:|
|none|0, t>=1|t|0|2|2|
|none|1<=q<2|t+q/sqrt2|0|2−q²/2|3|
|simple at8|0<=q<=sqrt2|t+1+q/sqrt2|1|1−q²/2|3|
|doubled at8|q=0|t+sqrt2|sqrt2|0|2|
|simple at7|q=0|t+1|2|0|4|

For example, a doubled contact at8 or a simple one at7 already has
projection cost2 onto B9, forcing q=0. A simple contact at8 has cost1,
forcing q²/2<=1. The continuous q intervals in the table are relaxed
necessary domains containing all native possibilities. The first row
has t>0 by positive height and hence t>=1 by the native product rule.

The source rederives all nine profile identities and all20 cross
identities from the physical C Gram. For each upper bound it multiplies
maximum−g by the positive height squared and checks a polynomial with
nonnegative coefficients after shifting t to its lower endpoint and,
in the tail-only row, q to1+u. Thus the inequalities cover every t and
q in the stated infinite domains, not merely sampled native labels.

Both sides have g<=4. Equation(1) would therefore have right side at
least1−8/9=1/9>0. A direct path1 or2 must occur after all, and supplies
a distinct full unit. Neither the original near-terminal payment nor
another changed affine payment is needed. The two new weights may be0,
so C has55−3=52 after retaining the earlier endpoint zero.

## B: retain all profiles instead of stopping at another affine half

Use the canonical affine B10 diagram with spine0,...,8, doubled bond
7–8 and short leaves9,10 at0. Set a=8,b=9. The new zero orbit consists
of(8,9),(8,10); each formerly had weight1/2. The resistance isR=10.

Deleting either marked leaf gives B10 or D10. The complete four-profile
classification in the coupled theorem applies: isolated, end3, end4,
next3. It follows from the endpoint inverse diagonals1,2 and larger
remaining diagonals, together with the same cycle/exceptional-type
exclusion. We retain all four profiles, including those that formerly
supplied another affine half.

At the heavy end every profile has cost at most4. At the fork end,
end4 and next3 also have cost at most4. Isolated and end3 reach8 only
at t=1 and t=0 respectively; each such extreme has height1/2 and one
simple contact to the old selected fork leaf or its mate. Every
nonextreme isolated profile has t>=sqrt2 and cost at most4; every
nonextreme end3 has t>=1 and cost at most3. The exact source checks
all8 physical profiles and their complete native-ray polynomial bounds.

If the fork is nonextreme, the total cost is at most8 and(1) is
impossible. Hence it is extreme. Let x be its spectator vector and
let y,k belong to the heavy root. Equation(1) requires

    E=(2x−y/k)^T G_Z^(-1)(2x−y/k)=g_heavy−2<=2.        (2)

If x is nonzero at i, disjoint supports and Gram-metric Cauchy give
E>=4x_i²/2>=2. Equality is necessary, forcing x_i=1 and

    2x−y/k=G_Z e_i.

The nonpositive off-diagonal signs make i the only x-coordinate. If
i has a neighbor j, then y_j>0 and the actual c–j–i–d path3 gives a
retained full unit. Otherwise i is a whole A1 component, all other
x/y contacts vanish, and the heavy cost is exactly4. The only heavy
profiles with that cost are isolated t=1 and next3 t=0. No rank bound
on Z was used in this argument.

If x is empty, choose the actual corner containing the extreme root
and the far old heavy mark. A contact to the omitted old selected leaf
leaves a whole A1 against a separate whole B10 at its doubled terminal.
The upgraded unit pays the entire face. A contact to the old mate
instead gives a whole finite B11 component and an unchanged unit.
The heavy root is absent from this corner, so arbitrary contacts y
cannot extend either selected component.

In the remaining equality case x is a single contact to a whole A1,
y=0. If the extreme root attaches to the old mate, its selected finite
component is whole B12, again paid fully. Otherwise its corner has
whole A2/B10 and contributes the existing half. A second distinct
corner supplies another existing half:

* heavy isolated t=1 leaves whole A1 against whole D10 at a short leaf;
* heavy next3 t=0 gives the whole affine D10 opposite-fork pair, which
  retains1/2.

The source constructs all four heavy-profile/old-fork-attachment
combinations, verifies every strict white/edge principal, and identifies
the two actual corner components, roles and positive affine kernel.
No omitted spectator can extend them: Z is orthogonal to F, and the
proved equality leaves exactly the one whole A1 and no other new contact.
These two halves are distinct from the original B corner. Their sum is1,
so the original B half is unnecessary. Both B weights may become0,
leaving55−2=53.

## Coupon cost and simultaneous priorities

The only finite change is whole A1/B10 doubled-terminal weight1/2→1.
There is one marked choice in A1 and one eligible terminal in B10,
so the actual increment is exactly1/2 per component pair. It must be
written into any subsequent actual-cost rule.

The inherited finite representative maps B10 to D10 and A1 to itself.
The old actual A1/B10 doubled cost is1; its new doubled cost is2.
The unchanged A1/D10 representative doubled allowance is also2, for
the two short-fork-leaf halves. The source checks these exact values
against the frozen current coupon function. Every other type pair is
unchanged. Thus every existing simultaneous representative upper
comparison survives this upgrade; no residual optimization is needed
for that pointwise statement. A subsequent combined forest certificate
must still record the actual increment and all other new local changes.

For simultaneous quadrilateral coverage, first handle any actual
new C near-terminal corner by its independent path1/2 unit. The earlier
C endpoint-zero theorem has the same independent output. Next handle
any new B zero corner by the argument above: its outputs are finite
B11/B12 units, path3 units, the upgraded A1/B10 unit, or two distinct
retained A2/B10 and A1/D10 or affineD10 halves. None is a changed C or
B pair, a long finite rank15–17 payment, or a changed affineD16 pair.

If neither new C nor new B corner exists, weights on all actual ideal
corners agree with the previous rule, and finite weights have only
increased. The complete old whole-face proof therefore applies
pointwise. This handles the old coupled D11 proof's occasional B-half
output without recursion: an actual such B corner is in the earlier
independent branch. The new order16 D105 theorem uses a whole D5/D11
finite half and retained short8; if that route reaches a B zero, its
independent proof here supplies the complete payment.

All new changed pairs have distance9. Finite7/ideal5 triangle payments,
protected ideal6/7 roles, whole affine order<=9 outputs and the cycle11
triangle exclusion remain unchanged. The whole finite outputs above
have rank<=12, or are separate-component coupons, and therefore avoid
all deep finite8 zeros and recent finite long-pair reductions. This
proves one simultaneous local rule with triangle minimum2 and
quadrilateral minimum1.

## Verification scope

Final source-bound generation `287a9e` and fresh saved replay `b625d5`
both exited0/PASS. A separate metadata audit passed all22 direct
input/self bindings, all6 report links and source syntax. The parent
and a second agent each completed a full source/report mathematical
read, including the inherited physical diagram/profile and marked
affine-type helpers. The parent also reread the old C endpoint theorem.
Both checked the complete B9+A1 contact domains, C strict gap, all B
profiles and native gaps, metric-Cauchy equality and whole-component
sealing, distinct corner alternatives, explicit coupon increment,
noncircular actual-corner priorities, and final full weighted-table and
native-involution guards. Neither reviewer repeated profile, determinant,
literal-corner, scalar or historical census arithmetic. No concrete gap
was found in the stated local theorem.

The source performs only the exact small profile,
determinant and literal-corner computations described above. Published
finite/affine classification, actual factor-fiber geometry and the
completed whole-face rule remain explicit inherited inputs. No native
catalogue, arbitrary-support census, ideal convolution or finite DP
has been rerun or claimed here.
Source, certificate and report are final/frozen.

    python3 scripts/fresh_order11_bc_zero_transfer_20261018.py
    python3 scripts/fresh_order11_bc_zero_transfer_20261018.py --verify
