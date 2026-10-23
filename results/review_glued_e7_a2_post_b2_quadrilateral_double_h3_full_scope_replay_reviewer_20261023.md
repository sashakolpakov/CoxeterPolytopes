# Full outside-author review: glued post-B2 return and double-H3 edge

23 October 2026. **PASS for the complete written scope stated below;
no source correction requested.** The reviewer read the entire317-line
candidate, its five literal inputs at their required interfaces, and
manually rederived the new geometry, original-module congruences, finite
profile classification, closed strips and symmetry diagnostics. This is
a mathematical source review, not an execution or implementation audit.
No worker/import, enumeration, symbolic calculation process, test,
compiler, Lean, TeX, plan, Sphinx or Git operation ran. Metadata reads
and hashes do not supply mathematical replay credit.

## 1. Exact source, direct pins and contribution disclosure

Reviewed source:

`results/compact_h9_glued_e7_a2_post_b2_quadrilateral_and_double_h3_edge_candidate_20261023.md`

SHA256 `151a4185d99fac4a5f27ae68301387e9ba1ff09d1eb23c03f1471a216c2970ce`,
14,062 bytes/317 lines. No source bytes were changed by this review.

All five literal pins match current bytes:

| Input | SHA256 | Bytes | Exact inherited boundary checked |
|---|---|---:|---|
| results/compact_h9_glued_e7_a2_ramified_five_module_and_first_edge_candidate_20261022.md | 33281068c43fda60f3356dc2cfcc226fcc72aae03d818ce49598e14382deb17c |14506| ORIGINAL E7-plus-ternary module, admissibility, primitive norm pool and full arrangement |
| results/compact_h9_glued_e7_a2_d4_leaf_actual_norm_four_b2_edge_candidate_20261023.md | 35f39b3556beab980cb506f9df557e55688b0d8bc02a14986d4311ed6f7c4d09 |12268| Actual primitive B, z, lambda and FULL H3+A3+B2+[n] starting germ |
| results/review_glued_e7_a2_d4_leaf_actual_norm_four_b2_edge_root_acceptance_20261023.md | b451da593fd6a42010c2cdec48140db1a090a0a42d918a895b06de80047ea795 |2690| Written acceptance of that exact B theorem, including its nonsimple-G tie resolution |
| results/compact_h9_glued_e7_a2_h3_d4_complete_right_hexagon_face_candidate_20261023.md | f54f3909f6fcab280eb2a7f36aa9631957cb8e79da1a36e27dfd6dfbdd2cac32 |13413| Q is an ACTUAL ambient facet and t is an actual finite vertex in the SAME chamber |
| results/review_glued_e7_a2_complete_right_hexagon_root_acceptance_20261023.md | 63a63ef6dd7938b03bcd6a39e0a2b8b8202f24488a72d78e4a6ae4050932229c |2749| Complete accepted hexagon interface, without ambient finite-volume or integral-period credit |

The reviewer did not author this new317-line construction. I contributed
earlier finite-edge/actual-germ and original-module methods, and reviewed
related supporting construction mathematics. These shared foundations
are disclosed; this is an independent review of the NEW composition and
new calculations, not a claim of independently originating or replaying
all inherited inputs. The author's inherited contributions and root's
separate H3-swap contribution remain disclosed in the candidate.

The accepted starting B is essential: G at the same old parameter is
NOT an actual facet. The candidate uses B throughout and never replaces
the actual B2 germ by the rejected orthogonal G germ. Likewise the face
normalizers in the hexagon input are not imported as ambient isometries.

## 2. Original module, normalized roots and starting data

I checked the use of the ORIGINAL direct O-basis f1,...,f7,e1,e2,g,
with O=Z[phi] and ternary Gram [[2,-1,1],[-1,2,0],[1,0,-2phi]].
Half-coefficients in fundamental weights are not extra basis vectors.
The denominator-three gluing in the orthogonalized description does
not permit arbitrary halves in this original basis.

The inherited even primitive norm alphabet is exactly r²=2d,
d=1,2,b,2b, b=2+phi, with original reflection condition r.L/d in O.
Its derivation uses primitivity and the inverse original Gram, not a
prime-only norm guess. All four classes are retained in the new domain.
The nonphysical positive-definite embedding, local finiteness of the
full integral-root arrangement, finite spherical classification and
finite-vertex simplicity are explicit inherited structural premises.
None requires an assumed finite-volume ambient chamber here.

The displayed W,theta,omega5,T0,F,lambda,B,z agree with the pinned
actual B theorem. The starting nine-wall Gram has rank9 and positive
determinant64b*phi^-2, while z²=-phi^-7/4. In particular the positive
line calculations below start from a FULL actual germ, not a root
list conjectured from known wall signs.

## 3. Q is uniquely first on the released-f7 edge

For eta=f7-B/2 I rederived eta²=1, eta.f7=1, eta.B=0 and all other
retained orthogonalities. The products z.Q=-phi^-2, eta.Q=phi³ give
the finite closed upper s=phi^-5 and zQ²=-phi^-13/4. The whole segment
is future timelike by continuity from z.

The exact identity Q=4phi^5*z+phi³*eta-B/2 holds in original K-space.
For an alleged distinct actual first wall r at0<s<=phi^-5, J=r.eta>0
and v=-r.B>=0, hence

    r.Q=J(phi³-4phi^5*s)+v/2>=J(phi³-4)+v/2>0.

This contradicts actual-wall acuteness because Q is already ACTUAL from
the accepted hexagon. The self case r=Q is explicitly exempted; neither
the norm of Q nor its upper bound is incorrectly tested against itself.
Starting retained walls have J=0, and released-f7 has the wrong outgoing
sign. Thus this shadow argument covers every other norm class and
unknown wall without a new census. Local finiteness gives existence,
and the inequality rules out every distinct earlier or tied facet.

The endpoint germ is H3+A3+B2(B,Q)+[n]. Its B2 product is B.Q=-2,
and its determinant is64b*phi^-2, as claimed.

## 4. The complete finite return quadrilateral and its limited symmetry

I expanded the return identity

    zQ+(phi^-6/4)(-B-Q)=(phi^-4/2)*t.

The direction -B-Q retains Q and has negative B product. Both endpoint
rays belong to the SAME convex chamber, so their segment stays in it.
An unlisted wall cannot cut its relative interior: a linear functional
nonpositive at both endpoints cannot vanish only in the middle. An
additional wall vanishing on the whole segment would contradict the
full finite endpoint germs. Thus this is the complete connecting edge.

The four established edges close the face retaining H3+A3+[n].
Convexity gives the whole quadrilateral, rather than merely four edges
of a possibly larger polygon. The projected-f4 normal is -lambda of
square1; the others are -f7,B,Q with squares2,4,2. Their consecutive
products give angles pi/2,pi/4,pi/4,pi/2. The Gauss--Bonnet area is
therefore pi/2. No arbitrary Coxeter-face inheritance is assumed.

For K0=phi^-1*T0+phi*f7, K0²=4 and K0.f3=1. Its metric reflection
fixes B,H3,A3,n and maps z to phi³*zQ, so it matches the indicated
future rays. It is not ORIGINAL integral: f3-K0/2 has nonintegral
original E7 coordinates. The only other oriented full-germ matching
adds the A3 endpoint swap, the reflection in F+f6 of square4.
Their product sends f3 to f3-[K0+phi(F+f6)]/2; the original g
coefficient is -(2phi²+phi^4)/2=-(5phi+4)/2, not in O.
I checked that H3 is rigid, A3 has just this endpoint swap, and the
B2 norm2/4 roots cannot be exchanged by a K-linear integral isometry
(that would require sqrt2 in K). The isolated norm2b component is
distinct. Thus these matching germs supply no original chamber period.

## 5. V is an original finite upper, not an acuteness premise

For mu=(3f6+2f5-F)/4 I rederived mu²=3/4, mu.f6=1 and its zero
products with the other eight retained walls and z. The vector
V=phi(g+omega5) has norm2, unit original g coefficient, and the
displayed complete original covectors. In particular it is primitive
and reflects the ORIGINAL module integrally.

Direct products give z.V=-phi^-4/2, mu.V=phi/2, hence closed upper
tau=phi^-5. At w=z+phi^-5*mu, I checked

    w²=-phi^-11/2, w.g=-phi,
    w.f3=phi^-3/2, w.f4=phi^-4, w.f6=phi^-5.

The old-f3/f4 functions stay strictly positive over the whole closed
interval, and released-f6 is strict for positive tau. They cannot be
self-wall exceptions to the first-wall acuteness signs h,l,p>=0.
The proof correctly does not assume V actual at this stage: an integral
upper plus the full initial germ and local finiteness supplies some
actual first facet, which must then lie in the complete domain below.

## 6. Complete finite profiles and the original half-class

A first wall touching retained H3 could only extend H4 by native3 at
-f2, with d1 and no other retained contact. Its E7 expression has
l+p in2O but R orthogonality gives2x=3+4h+3l+p. With n orthogonality,
x must be in O, contradicting the odd constant. This includes both
released directions, so no hidden H4 branch is lost.

For the remaining profiles I independently rederived all equations(7),
including the ORIGINAL conditions l+p+k in2O,
2x-(phi-1)j in2O and2y-j in2O. Their intersection gives
j=phi*k and p+t=phi²*k modulo2O. The exact B equation is
-v=-2phi^-1*h+l+phi^-2*k-t. These are original membership constraints,
not an assumed integral individual Cartan-entry normalization.

For retained n, all nonzero finite options are precisely j1 at d1,
j1/phi at db and j2 at d2b; d2 has none. A native10 cannot extend
to finite rank3. The possible unit-root contact of a db wall would
itself be native10; a d2b wall has no such finite contact. Thus these
rows have no other retained contact. The j1/phi rows violate the mod2
equations. In d2b,j2, reflection divisibility forces h,l,p in2bO and
x,y in bO, giving r=b*s-n with ORIGINAL s, s.n=0. Its norm would
require b*(s²/2)=1 in this even module, impossible. Therefore j=0
for every surviving first wall, without deleting a foreign n contact
by assumption.

The B2 inverse-Gram cost k²+kv+v²/2 is correct. The full finite list
leaves short-side3 (d1,k1) or long-side3 (d2,v2). All other listed
bonds are null/indefinite or require square roots outside K. The
short-side3 and its possible A2 bridges fail k parity. The long-side3
cannot join A2: projection cost2+8/3>4. The A2-alone possibilities
are d1 native3/5 at either end or d2 native4 at either end. Double
contacts create a nonpositive circuit, and native10 cannot extend A2.

This gives exactly ELEVEN necessary profiles: four orthogonal norm
classes, two d1/native3, two d1/native5, two d2/native4, and one
d2/v2. Both A2 orientations and all bridges were checked separately;
no presumed original-module diagram symmetry is used. The residual
C is2d-(2/3)(q²+qt+t²)-v²/2 and has precisely the displayed values.

## 7. Complete closed-strip proof and unique V survivor

After k=j=0, l=p=t modulo2O makes the displayed time aF even; thus
Y=aF/(2phi³), Dold, X and D really lie in O. I rederived

    r.z=-phi^-2 X/2, r.mu=D/2,
    D²-3phi³ X²=3C,
    X>0, D>0, D>=phi³ X.

The last condition includes the CLOSED upper. The physical coefficient
identity phi^6-3phi³=2phi² and the conjugate norm equation give exactly
the two strips(9). No unbounded time/Pell domain is being substituted.

Orthogonal b/2b rows have whole original b-content after covector and
evenness conditions, hence are imprimitive. Orthogonal d1/2 and the
pure norm4/v2 row reduce modulo3 to X²=phi^5 or phi in F9. Phi has
order8 there, so both are nonsquares. This excludes all isolated norm
classes and the B2 long-side row, including both signs of a prospective
field square root.

For C4/3 I checked the closed physical/conjugate bounds, the exact
coefficient reduction B=-1,0,1 in X=A+Bphi, and the sole strip value
X=phi-1. Its D²=4+3phi has rational norm19, not a rational square;
both native3 and both norm4/native4 orientations are excluded.

For C=2phi^-2/3 the only strip value is X=2-phi, including its closed
physical equality. Then D²=phi² and the positive collision sign selects
D=phi. The t=phi orientation forces p=phi/3, nonintegral. The q=phi
orientation forces p=0,l=2L,h=phi L. Its exact starting product is
-phi^-5 L-phi^-4/2, so X=phi^-2 forces L=0. The remaining original
normal is exactly V. Negative D cannot cut the outgoing edge; X=0
cannot be a new opposite from the FULL starting germ. Thus no closed
upper case or tie was silently omitted.

Existence of an actual first wall plus this exhaustive necessary-domain
proof forces V actual and uniquely first. Unknown-wall feasibility
is derived from firstness, rather than inferred merely from the named
wall signs or from V's membership in the integral root arrangement.

## 8. Full double-H3 corner, swap diagnostic and exact conclusion

The new nine actual normals have full positive Gram
H3(R,-f1,-f2)+H3(V,-f5,F)+B2(-f7,B)+[n]. The new H3 bonds are5,3;
released-f6 is strict and is not mistakenly included in that H3.
The determinant32b*phi^-4 and future finite ray w check. Simplicity
therefore proves this is the FULL ambient germ, not just a subdiagram.

For the proposed H3 exchange, D1,D2,D3 are the old-minus-new oriented
normal differences, with Gram2G. I checked the projection formula
J(x)=x-D G^-1(D.x), D.f3=(0,0,1-phi), and equation(11).
The resulting g coefficient is exactly -(5phi+3)/2, not in O.
Since both H3 diagrams are rigid and the other norm classes/components
cannot be exchanged, this excludes the only nontrivial permutation of
this oriented full germ in the ORIGINAL module. This statement concerns
oriented chamber-germ symmetries; it does not deny the usual finite
reflection group fixing the hyperbolic vertex while changing sectors.

**Accepted review scope:** two complete actual edges from z, the actual
finite return quadrilateral of area pi/2, the unique new primitive V
facet, its full double-H3 corner, and the stated NONintegrality of the
matching-germ maps. No mathematical defect was found in that scope.
This establishes no ambient H9 finite volume, chamber period,
reflectivity/nonreflectivity, new compact example, field/form-class
exclusion or dimension improvement. Root acceptance remains a separate
gate; no execution or formal-verification status is promoted here.
