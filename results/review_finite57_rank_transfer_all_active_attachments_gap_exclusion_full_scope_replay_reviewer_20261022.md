# Full mathematical review: the stronger heavy-transfer third-corner bound

22 October2026. **PASS at the exact restricted scope below.** This is
a complete manual mathematical source review of the root-authored

`results/finite57_rank_transfer_all_active_attachments_gap_exclusion_candidate_20261022.md`

SHA-256 `c6430be1f294e6b784a514675af1764c16569baaefa2d952b5cafb1a60491d00`.

I read the entire proof and the necessary pinned inputs, reconstructed
the finite extension lists, all six active orientations, every new Schur
calculation and all charge maxima. No proof worker, import, enumeration,
test, Lean or TeX job ran. No source, guide, paper or bound was edited.

## Dependency and contribution boundary

I authored the inherited heavy-transfer source8abe and contributed to
earlier Context662 financing and source audits. This is an outside-author
review of the NEW root-authored exclusions and their complete composition;
it is not presented as an independent referee report of my own foundation.
The current calculation does not spend any inherited selected-recipient
surplus or create an allocation. Input bytes were streamed and matched:

| Source | SHA-256 | Scope read |
|---|---|---|
| results/finite57_heavy_rank_transfer_quadrilateral_local_discharge_20261021.md | 8abe292196a6792ed350c85c0c1e66ce397f312af2a025c2502a287a601f3f53 | Actual full white support, all d extensions, previous three heavier-row exclusions, exact C0 convention |
| manuscripts/finite_volume57.tex | c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0 | Actual-face positive projection/opposite-core lemma; current finite same/cross/context rules |
| results/finite57_d16_full_orbit_priority_with_retained_context662_credits_20261021.md | 0b09e32f10721337ddff468dc521aa1a31178d92d7b1e619747ef9113fdd01ae | Held order11 increments and selected-credit versus actual-charge distinction |
| results/finite57_homogeneous_large_d_quadrilateral_schur_obstruction_20261021.md | 3ef4130fc30f557425fdb06144342054e019c6eedda3e54a197c532b158b10ac | Full-white Schur/rank method and actual norm-two native gap |

The latter two are inherited dependencies of8abe. The new proof directly
pins8abe and c300. Acceptance labels were not used in place of their
mathematical interfaces.

## 1. Exact geometric statement and exhaustive finite alternatives

The theorem assumes an ACTUAL quadrilateral in H57, full white forest
W=A13+2A1+D17+D23, cyclic sides a,b,c,d, the specified long-end/center
roles at(a,b), and c's sole white simple3 contact at the inert D23 long
leaf. Both d corners are finite. All roots have norm2. These are essential
hypotheses, not an asserted realization of a formal heavy forest.

I checked the finite extension classification. At(c,d), W+c contains
A13+2A1+D17+D24. A finite attachment to a large D factor can only
extend its long arm by simple3. It cannot join both large branched
factors, add a second fork, or introduce a doubled edge into a large
branched finite type. Its only additional contacts are simple3 to ONE
isolated leaf or to ONE endpoint of A13; doing both creates another
fork. This gives exactly the six listed forests. In particular an
attachment to c is an extension of the D24 long arm and is included,
not discarded by setting c.d=0 prematurely.

If neither large factor is joined, the remaining finite list is exactly
A13+3A1, A13+T2+A1, A13+T3, A14+2A1, B14+2A1, D14+2A1,
A15+A1, B15+A1, D16. The small types are whole irreducible factors
and include the finite dihedral and H3 possibilities. An attachment
to A13 is endpoint3/4 or next-to-end3: two contacts form a cycle,
and no rank14 exceptional finite diagram supplies another case.
Joining an isolated leaf extends the endpoint path; joining both makes
the indicated D16. The inherited heavier three forms have only their
far orientations compatible with the finite(a,d) corner; the full-white
native-gap exclusions in8abe cover them. Their use here is legitimate.

Thus the NEW remaining work is precisely the six orientations of
A14+2A1, B14+2A1, D14+2A1. Here d is orthogonal to c because it
does not join the D24 block. This orthogonality is conditional on this
branch of the exhaustive list.

## 2. Orientation check and exact rank equation

The native inverse for A13 has diagonal j(14-j)/14 and endpoints
inverse product1/14. Each isolated white has inverse diagonal1/2;
the D23 long leaf has inverse diagonal1. Eliminating ALL55 whites gives

    S_aa=15/14, S_bb=1/14, S_cc=1,
    S_ab=-1/14, S_ac=-g, S_bc=S_cd=0.

These include both old isolated leaves in b's norm and the entire
inert D23 in c's norm. There are59 actual roots in a58-dimensional
Lorentz space, with the55 whites positive and independent. Therefore
the4-by-4 Schur determinant is0. Eliminating the positive c coordinate
leaves a's diagonal15/14-g².

For all six orientations, any nonzero RAW a.d would create a cycle
through the connected A13 path at the finite corner(a,d). Therefore
raw a.d=0 before elimination, not after it. The two immediately
nonfinite near orientations are correctly removed:

* endpoint4 at q1 has the principal star at q1 with edge magnitudes
  1,sqrt2,1. Its central Schur square is2-(1+2+1)/2=0, so it
  cannot be a principal subset of the finite corner;
* next-to-end3 at q2 has arms2,1,11 and rank15, outside the finite
  simply-laced A/D/E list.

The near endpoint3 at q1 DOES instead give a finite D15 at(a,d).
The new source retains it for a separate proof. The orientations at
q13/q12 survive finite classification and are handled by Schur.

For actual distinct norm-two wall normals, native finite labels give
g² in{0,1} union[2,infinity) and h=0 or h>=1; disjoint/parallel
cases also lie in this domain. This uses the actual positive-principal
face correspondence, not Coxeter angles assigned to a projected face.

## 3. Every far Schur exclusion rederived

Write B0 for the b,d block and u=-1/14 for S_ab. With v=S_ad,
detB0 nonzero, the rank equation is exactly

    g²=15/14-(u,v)B0^-1(u,v)^t.

Here is the complete recomputation:

| d attachment | S_dd | v | -S_bd-h | -detB0=Q | Result |
|---|---:|---:|---:|---|---|
| q13 simple3 |15/14|-1/14|13/14|h²+13h/7+11/14|g²=15/14+(3+2h)/(196Q)|
| q13 simple4 |1/7|-sqrt2/14|13sqrt2/14|h²+13sqrt2 h/7+12/7|g²=15/14+(2+sqrt2 h)/(98Q)|
| q12 simple3 |2/7|-1/7|6/7|h²+12h/7+5/7|g²=15/14+(h+1)/(49Q)|

Every Q is strictly positive on ALL h>=0. Multiplying Q respectively
by42/11,7/6,7/5 and subtracting the relevant numerator gives

    (42/11)h²+(56/11)h,
    (7/6)h²+(7sqrt2/6)h,
    (7/5)h²+(7/5)h.

Thus the three exact intervals are respectively

    15/14<g²<=12/11,
    15/14<g²<=13/12,
    15/14<g²<=11/10.

All are strictly inside(1,2), contradicting the ACTUAL native gap.
No finite, parallel or dotted value of the opposite h is omitted, and
no unmentioned numerical extremum or square-root approximation is used.

## 4. Near endpoint3, including the essential h=0 boundary

For q1 simple3 the correct entries are

    S_dd=15/14, v=-13/14,
    B0=[[1/14,-h-1/14],[-h-1/14,15/14]],
    detB0=(1-2h-14h²)/14.

At h=0 both diagonal and determinant tests are positive. Consequently
the ORIGINAL W+b+d principal set is positive definite. The c300 positive
projection constructs a common finite point feasible against every actual
wall and lying in the original face defined by W. This would make the
opposite closed sides b,d meet, forbidden for this actual nondegenerate
quadrilateral. This removes h=0; a Schur/native-gap argument alone would
not remove it. The remaining native domain is h>=1, where B0 is invertible.

The inverse calculation gives exactly

    g²=h(15h+4)/(14h²+2h-1).

The positive denominator subtracted from its numerator is(h+1)²;
twice the denominator minus the numerator is13h²-2>0. Therefore
1<g²<2 throughout h>=1. This is the required contradiction and
completes all six orientations.

## 5. Complete charge maximum and scope

I recomputed the six large-D values from the actual c300 same/cross/
presence rule. The whole D costs used are D18=112,D19=120,D24=150,
D25=155,D26=160,D31=185,D38=220,D17=113.5. A13 has ordinary77,
increased to78 when a whole A1 remains; its A1 cross total is137/64.
No A1 has a same-pair cost. These give exactly

    344.28125,350.78125,350.140625,353.640625,335,333.5.

For the untouched path the active values are84.421875,85.859375,
81.5. The middle maximum is at T2=A2: same cost1, cross151/32
with A13 and the A1 cross137/64. Thus its full value is
263.5+85.859375=349.359375. Other rank2 types have no A2-specific
endpoint upgrade. For T3, the full3-node same cost3 and broad A13
center cross1.5 give81.5, with no whole A1 guard.

Every contextual guard relevant to these exact forests is included;
the held order11 increments all vanish because no B8/B10/D10 occurs.
No late selected Context662 credit is subtracted from an unselected
actual finite vertex. After the old three heavier forms and the new
three active forms have been excluded, the exact maximum is therefore

    C0(c,d)<=22633/64=353.640625.

The illustrative mean(372+22633/64)/2=46441/128=362.8203125
and the dimension57 coefficient364=(56*52)/8 are correct. This is
only a numerical opportunity for an allocation. The proof supplies no
global pairing, injection, recipient multiplicity control, treatment
of different c roles, ideal endpoints or nonquadrilateral faces.

**Final verdict:** full manual PASS of the stronger third-corner bound
under EXACTLY the stated actual quadrilateral/whole-white/c-contact and
finite-endpoint hypotheses. It is pending root acceptance as a new
supporting theorem. It neither deletes a whole coupon nor improves the
global charge372 or the finite-volume ceiling57, and it asserts no
example or runtime/formal verification. No correction is requested.
