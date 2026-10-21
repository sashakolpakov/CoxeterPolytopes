# A retained-wall conic and a new six-cycle cavity at wall66

This bounded numerical experiment keeps **all63 old walls**, including wall8,
and changes the wall66 region by adding a cap on the single old59/60 edge.
It finds78 caps passing the native and common-center tests. Three distinct
low-label caps expose a new six-cycle and admit274 further single caps, but
all6,119 compatible pairs are dotted. Every accepted companion single or pair
still has **ten open ends**. No completed disk, exact new normal, finite-volume
chamber or new-dimensional example is claimed.

The [source](../scripts/fresh_wall66_single_edge_conic_20261018.py),
[saved certificate](fresh_wall66_single_edge_conic_20261018.json), and
[independent saved-record checker](../scripts/check_wall66_single_edge_conic_saved_20261018.py)
are new. The exact
[63-wall control](fresh_nonarithmetic_release49_62_replacement_20261017.md)
and its [complete S6 boundary](fresh_nonarithmetic_release49_62_boundary_20261017.md)
remain unchanged:199 complete ambient finite vertices, ten target-face open
ends and a certified escape, hence infinite volume. This experiment continues
the search for a closed chamber in H9; its S6 cone is only a three-dimensional
face of that ambient chamber.

## Changed duty and deletion ledger

Earlier [wall66 release](fresh_release66_three_cap_disk_20261017.md) and
[joint58/66 release](fresh_release58_66_seven_boundary_20261017.md) attempts
required caps on two adjacent old boundary edges. Their empty banks and later
shifted-anchor obstructions are not rerun here. The new initial capK is
required to form just the finite corner59/60/K. It has a continuous degree
of freedom before any third old-wall native contact is imposed.

There is no wall deletion: the debit is exactly zero old finite supports,
zero compact faces and zero supports involving wall8. This is a deletion
statement about unchanged old constraints, not an exact incidence certificate
for numerically discovered additions. Every new normal is screened against
all63 old ambient walls, through the inherited S6 profile and native-offset
representation. No old contact, in particular58–59, is released.

The [exact old58/59 construction](fresh_nonarithmetic_pi8_quad_completion_20261017.md)
prescribes their native label9. This is not a rounded label inferred from
the scout. A norm-two finite triple containing that bond cannot have another
nonright bond: with theta9=2cos(pi/9), its determinant is at most
8−2(theta9²+1)<0. Thus any finite58/59/K triple would requireK right to both.
For profile0, a cap right to58,59,60 cannot have norm2, since those three
intrinsic normals span a positive-definite three-space in Lorentz4.
The chosen seeds avoid imposing the58/59/K corner. They create the other
two endpoints of the old60 face instead.

## Two-sheet parametrization and native events

Work in the rest frame of the inherited common center. LetU be the
positive-definite span of the intrinsic59 and60 normals. For a fixed A2
profile and prescribed native contacts, letp be the unique vector inU
with the required intrinsic products after subtracting native offsets.
The saved projected norm gives d=q−p²>0. In the Lorentzian complement ofU,
choose a future unit timelike vectorh and a unit spacelike vectore, with
h·e=0. Both sheets of the norm equation are

    K(t,s) = p + sqrt(d)/2 [s(t+1/t)e + (t−1/t)h],
    t>0, s in {−1,+1}.

Every native product is a*t+c+b/t, including the fixed profile offset.
Acuteness, strict dotted contact and the common-center condition therefore
reduce to quadratic sign inequalities after multiplication by t>0.
A third-wall native labelm is the quadratic equation

    a*t² + (c+2cos(pi/m))*t + b = 0.

The inherited finite-corner classifier gives the symbolic59/60 recipes
`(profile,m59,m60)=(0,2,2),(0,3,2),(0,4,2),(0,5,2),(0,H,2),(1,2,2)`.
HereH denotes an arbitrary label at least6. The present experiment samples
H=6,8,12,24, giving nine recipes and18 sheets. Third-wall events use
labels2,3,4,5,6,8,12,24. All57 projected old normals are tested; their
six S6 products are encoded by the seven inherited profiles.

The code intersects numerical quadratic sign intervals and solves the
native-event quadratics directly; it does not sample a grid in t. It also
tests one interior point from each surviving strictly dotted interval.
This does **not** certify every real parameter: interval endpoints,
near-degenerate coefficients, near-tangent cases and unsampled labels are
not enclosed by exact arithmetic. A recorded unresolved list is fail-safe
bookkeeping for encountered numerical conditions, not a proof that none
exists elsewhere.

The18 sheets give7,920 event equations and7,925 positive-parameter tests.
Of these,7,847 fail a native or center test. The78 distinct accepted caps
comprise nine strictly dotted interval representatives and69 events with
an additional finite old58 contact. No other old-wall native event survives
the stated search. Negative sheets have empty computed acuteness domains;
positive sheets have computed dotted tails. The caps have10,11 or13 open
ends, never fewer than10. Every acceptance rechecks the norm/contact
residual, all retained native products, strict center and the complete
numerical cone against an independent positive-definite-support ledger.
Finite labels use tolerance min(2e−7, one hundredth of adjacent cosine
spacing); a near-parallel product is not rounded to a large finite label.

## Three selected seeds and the six-cycle disk

The bounded companion test chooses the following profile0 seeds, all right
to old60. These are a selected domain, not all69 native events.

|Saved capK index|Labels to58,59,60|X bank|Y bank|Z bank|
|---|---|---:|---:|---:|
|1|(3,2,2)|6|41|36|
|9|(2,3,2)|0|119|36|
|10|(3,3,2)|0|0|36|

For each seed the complete numerical cone has28 finite and6 exterior rays,
with ten open ends. Its old60 face has exactly the three finite supports
`(58,59,60),(58,60,K),(59,60,K)` and no other rays. The six one-open edges
in the wall66 region form the verified cycle

    30–59–K–58–66–33–30.

The four independent B-region open edges remain7/47,7/30,47/52,30/52.
This separation is checked from the full saved edge ledger, not inferred
from a proposed drawing.

The companion frames areX=(30,59,K),Y=(K,58,66),Z=(66,33,30). Each cap owns
its two consecutive boundary edges. Three junction trianglesKXY,66YZ,30ZX
and the centralXYZ triangle complete the proposed disk. It hasV9/E18/F10;
every interior edge occurs twice, every boundary edge once, and all links
are connected intervals or cycles. The three cap links are pentagons.
The outer anchors of a frame are not falsely required to meet: the finite
classifier receives only the two actual boundary-edge duties.

The nine banks test486 norm equations and972 branches across all seven
profiles, low finite labels and the same four high samples. They accept274
single caps. All have30 finite rays,6 exterior rays and ten open ends.
The6,222 cross-bank pair tests give103 nonnative products and6,119 dotted
products. Every compatible pair is scored even though its dotted contact
cannot supply the required finite junction; all have32 finite rays,
6 exterior rays and ten open ends. There are no finite junctions and no
three-cap candidates. Full triple/native/center, ten-corner and whole-cap
face checks remain in the source but are not reached.

For a concrete boundary example, seed1 has exterior supports
`(7,30,47),(30,33,66),(30,47,52),(30,59,K),(30,66,K),(58,66,K)`.
After its first X/Y compatible pair, the exterior supports are
`(7,30,47),(30,33,66),(30,47,52),(30,66,X),(66,X,Y),(K,X,Y)`.
The new finite supports then are
`(58,60,K),(59,60,K),(30,59,X),(59,K,X),(58,66,Y),(58,K,Y)`.
Thus more finite vertices move the exterior region without closing it.
All18 distinct complete ledgers, including the other seeds and pairs,
are saved rather than replaced by only their counts.

The result supports changing the mutual companion geometry or duties in
any future attempt. It does not exclude the continuous single-edge family,
other native-event labels, other seed choices, or different cavity disks.
No further bank was run after these three authorized seeds.

## Verification and limits

Generation95375 terminated exit0/PASS (terminal87e83b), and fresh saved
replay11616 terminated exit0/PASS (terminald8fb39). The separate
standard-library checker terminated8dcf2f/PASS and final31406e/PASS:12 bindings,18 saved
ledgers,78 conic caps,274 companion caps and all6,222 pair records/accounting.
It rechecks saved conic coefficients/norm residuals and native values,
the exact inherited label9 declaration, zero deletion, disk links and
saved cone/PD support equality. It imports neither the producer nor a
numerical package and does not reconstruct physical normals, roots or cones.
All12 source bindings, both Python ASTs and eight local Markdown links
also pass the final metadata check. No owned job remains active.

An independent peer checked the literal six-cycle disk and its links;
that review is topology only. The inherited physical reconstruction,
spherical classifier and cone routines remain stated dependencies.
The root reviewer read the complete producer, saved checker and report,
and passed the mathematical/source scope: conic parametrization, native
event and sample limits, inherited label9 obstruction, zero-deletion
meaning, six-cycle duties and full pair-ledger accounting. That review
did not rerun numerical normals, roots or cones and does not independently
certify the inherited classifier. The packet is final and frozen.
All new results here are numerical discovery and saved-record checks.
The exact63/199/10 infinite-volume control remains authoritative, and
there is no full ambient finite census or volume certificate for an addition.
