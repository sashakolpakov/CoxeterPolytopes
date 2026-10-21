# A longer compact10 chain forces its first cap to be old44

On the frozen exact63-wall scaffold, delete30,33,44 and retain the other60
walls. Consider a first replacement J with finite corners8/10/J,
8/65/J and63/65/J, using the inherited seven possible white profiles.
Then every native J is the old wall44. There is no changed first cap
with these duties, at any finite Coxeter label. This is a necessary local
obstruction, conditional on the frozen exact scaffold and complete
seven-profile spherical classification. It is not a general chamber
nonexistence theorem or a finite-volume construction.

The [producer](../scripts/fresh_compact10_three_chain_gate_20261018.py),
[certificate](fresh_compact10_three_chain_gate_20261018.json) and
[independent interval checker](../scripts/check_compact10_three_chain_gate_saved_20261018.py)
preserve the [moving33 experiment](fresh_moving33_compact10_bridge_20261018.md).
The proposed longer wall10 boundary is7–11–8–J–R–S, retaining three
replacement sides rather than shortening the old boundary. The present
test stops before constructing S or R because J cannot change.

## Duties and exact deletion debit

The original compact target faces include8/10/44 and8/44/65 on wall8,
and44/63/65 on wall65. Replacing44 by J while retaining these incidences
therefore requires precisely the three finite corners above. There is no
assumed10/65/J corner. Later duties involving R and S have not been
tested here.

Deleting30,33,44 removes18 of the199 exact old finite-support
certificates,12 compact two-face certificates and one compact three-face
certificate. Two removed finite supports contain wall8. This is the
same exact certificate debit as the moving33 packet. It is not a count
of all faces of an accepted replacement chamber: no such replacement
exists in this family. The frozen original chamber still has ten target
open ends and remains infinite-volume.

## Original physical intervals and the norm line

Every old physical normal is reconstructed from its exact facet witness
and common-center product. All six white roots are projected out by
their positive-definite Gram matrix. Thus the projected ambient space
has dimension4 and signature(3,1). Original native equalities and the
seven white profiles are inherited from the exact control; interval
containment checks their physical interface and is not used to infer
an exact equality from a small floating-point value.

For a projected three-frame with Gram G and prescribed product vector
r, the frame signature is(2,1). Its orthogonal complement is a positive
line. If the proposed projected norm is q, write R=q−rᵀG⁻¹r. For a
retained gate wall g, write its projected frame coordinates as G⁻¹c and
its residual square as s=g²−cᵀG⁻¹c>0. Including the native white offsets,
the two possible products with g are

    L ± sqrt(R*s).

All displayed intervals use outward rational rounding. A continuous
domain is intersected with the necessary R>=0 condition, including the
zero-remainder boundary. No numerical normal is promoted to an exact
wall by this argument.

## All low recipes

The two initial finite duties8/10/J and8/65/J give46 low recipes over
all seven profiles and labels2,3,4,5. The producer also compares this
list to the complete7×4³ cube filtered by those two finite-corner tests.
Both norm signs are enclosed for each recipe.

The retained gates4 and63 exclude all92 signs:90 products are strictly
positive, and two lie strictly between−1 and0, which contains no native
norm-two Coxeter product. The residual cross product of gates4 and63
is strictly negative. Since their residuals lie on the same positive
line, their radical signs are opposite. This sign relation is checked
explicitly before using different gates to exclude the two branches;
independent sign choices are not silently identified.

## All five high-label families

The inherited complete finite classifier leaves exactly five symbolic
families. Four have profile0, J8 right, J10 label2,3,4 or5, and J65 label
at least6. The fifth has profile0, J8 label at least6, and J10,J65 right.
No high J10 bond is possible because old10 already has a simple contact
to a white root in the required finite corner. A high bond together
with another nonzero native bond makes that three-root principal Gram
non-positive-definite, including equality at labels6 and3.

For the fifth family, retain the entire closed magnitude interval
[sqrt(3),2], which contains every finite label at least6 and its limiting
value. In the frame8/10/65, the two J63 products are enclosed by

    [0.684709399, 2.048374031],
    [1.773943946, 3.137608579].

Both are strictly positive. These decimal bounds are outward displays
of the saved rational bounds, not rounded coefficients used as proof.

For the other four families, the additional finite corner63/65/J is
essential. Old63 and65 are right. If t>=sqrt(3) is the J65 magnitude
and x is the J63 magnitude, its principal determinant is
2(4−t²−x²). Strict positivity and the native alternative x=0 or x>=1
force x=0. Hence solve J first in frame8/10/63 with products
(0,−theta(m),0), m=2,3,4,5, and measure its product with65 afterward.
This is a finite native-first reduction of the whole high-label tail,
not a label cutoff.

The positive norm branch always has positive J65 product. For the other
branch, m=3 lies strictly between native labels12 and13; m=4 and5 lie
strictly between labels11 and12. The proof uses the exact theta12 value
(sqrt(6)+sqrt(2))/2 and elementary rational bounds

    theta11 < 2−(3/11)^2+(3/11)^4/12,
    theta13 > 2−(22/91)^2.

These follow from3<pi<22/7 and the alternating cosine bounds. Their
strict interval separations are saved. Thus these three branches are
nonnative regardless of which high label was proposed.

For m=2, all three frame products and all six white products vanish.
The one-dimensional positive complement is spanned by the exact old44
normal, which already has square2 and those exact right contacts.
Consequently the two normals are exactly old44 and its negative. The
strict common-center sign retains old44; its opposite also fails the
retained65 acute gate. This case is an original-wall control, not a
new algebraic normal or an approximate equality.

## Consequence and verification scope

The attempted longer chain cannot change its first wall while all three
listed duties and retained walls survive. A future construction must
change at least one of those incidences, the retained anchors, or the
allowed white-profile premise. Moving another later cap alone does not
evade this obstruction. This packet neither tests those alternatives
nor asserts compactness or finite volume from a target-face condition.

Final generation session60323/chunk8b515b and separate fresh replay
session59550/chunk824735 terminated exit0/PASS. The independent checker
passed in chunkfcf10d. The root agent completed a full producer, report
and helper mathematical/source read: physical projection, correlated
radical signs, low cube and inherited five-family scope, continuous norm
domain, native-first high65 reduction, rational angle bounds and exact
old44 complement/center argument passed. That read did not execute the
code or independently reconstruct the old scaffold, spherical classifier
or cones. Source, data, report and helper are final/frozen.

The independent checker imports only the standard library. It reconstructs the original physical interval Gram,
uses exact midpoint inversion with a Neumann-series error enclosure,
independently checks all92 low-sign exclusions, the continuous high8
box, four native-first high65 cases, sign correspondence and old-support
debit. It checks the saved46-recipe index domain but inherits its
spherical-classifier completeness. It does not construct normals,
replay cones or independently prove the exact old scaffold.

    python3 scripts/fresh_compact10_three_chain_gate_20261018.py --verify
    python3 scripts/check_compact10_three_chain_gate_saved_20261018.py
