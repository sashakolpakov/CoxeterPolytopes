# Independent review: the B19+A1 half-unit end refund

19 October 2026. Full source/report mathematical read, with independent
hand derivations below. No scientific job, matrix calculation program,
enumeration, optimization, Lean or build was run. The reviewed report is
unchanged at SHA256
`7d9b7ed96c0ddc56b225178b1920f30120b05218c9c91f41e433c901c7f20882`:
`results/arithmetic_b19_a1_spectator_half_end_refund_20261019.md`.

**Verdict: PASS as a written consumer theorem in the stated eight fields**,
relative to its explicit inherited arithmetic, finite-classification and
acute-face inputs. No missed continuation, sign, field or payment case was
found. This review performs no formal or computational verification and
does not promote a global dimension19 theorem.

## Full spectator continuation and Schur signs

At a selected vertex the complete white graph is B16+2A1, with isolated
b,h. In a newly unpaid single-source quad the other marked pairs have
original weight zero, so c·E=a·D=E·D=0. A contact E·b would give the
retained c−b−E distance-two unit. The actual finite diagram W+c already
contains B18, so every further contact into that large component must
extend its unique simple endpoint b. Since that contact was excluded,
E is either isolated or attached only to h with finite magnitude gamma.

Likewise D·b is excluded by a−b−D. The complete possibilities are H:
contact to h only, or T: simple contact to the simple endpoint S of B16,
possibly also a simple contact to h. In T the latter magnitude is0 or1/2;
any other bond would fail the large finite B-path classification. This
includes the possible spectator attachment and does not assume h stays
isolated at later corners. Simultaneous E,D contacts to h would give a
retained distance-two unit, so gamma*delta=0.

Eliminating W gives a_bar²=3/4, c_bar²=1/4 and a_bar·c_bar=−1/4.
The remote B16 endpoint inverse is2. Thus the T correction to c·D is
−1/2, while the E,D projection cross term vanishes. With

    t=1−gamma²,
    (u,Y)=(1−delta²,y) in H,
    (u,Y)=(1/2−delta²,y+1/2) in T,

the stated four-by-four Schur matrix is correct. Strict opposition gives
x²>3t/4 and Y²>u/4: a positive-definite opposite principal would project
to an actual intersection; a positive-null one would provide an ideal
point, contrary to compactness. Its determinant reduces exactly to

    (x²/t−3/4)*(Y²/u−1/4)=1/16.

This uses the actual complete white support and n+2 normals in Lorentz
dimension n+1, not an isolated three-node subdiagram.

## Arithmetic field step and complete scalar exclusions

The old20-root vertex span is positive definite. The new E wall is not
incident to that vertex, so its old-span projection has cost strictly
greater than1. Inverting the old a,c Schur block gives[2,2;2,6]; the h
block is orthogonal. Consequently q_E=2x²+gamma². At every nonphysical
embedding the independent old-plus-E block is positive definite, so
q_E'<1. These are determinant-ratio statements in the ground field;
they do not require individual unit-normal square roots to belong to K.

The squares x²,gamma² are cyclic quantities in K. If they generated a
proper subfield, a nonidentity embedding fixing that subfield would
preserve the physical q_E>1, contradicting strict conjugate positivity.
Thus K=Q(x²,gamma²). This step is essential to the exclusions.

For gamma!=0, delta=0 and the second determinant factor is at least1/4
in both H and T. Hence x²<=1−gamma²<=3/4, restricting x to labels3,4,5,6.
Combining this with q_E>1 leaves exactly the listed spectator labels:
(m,ell)=(3,5),(3,6),(4,3),(4,4),(5,3),(6,3). Their generated fields
are Q or Q(sqrt5), outside the eight fields. No unrestricted spectator
label is omitted by this finite reduction.

For gamma=0 and finite x, x²>3/4 and every nonphysical conjugate is in
(0,1/2). Field generation and degree<=5 give exactly the listed14 labels
above6. The eliminated conjugates use integers coprime to m; their ratios
a/m<1/4 supply the prohibited conjugate greater than1/2. The only
remaining labels are12 and18. Their determinant value A lies strictly
between1/2 and(3+sqrt3)/8<3/5. Each H/T subcase is then excluded:

* H,delta=0 puts y² in(1/2,3/5), a native gap.
* H,delta!=0 forces y=1/2 and delta² in(1/2,7/12), also a gap.
* T,delta=0 forces0<y<1/2.
* T,delta=1/2 forces(y+1/2)²<1/4.

For dotted x>1, the T determinant product is strictly too large. In H
it forces y=1/2, then delta=1/2, u=3/4 and x²=3/2. The rational value
q_E=3 contradicts every nonidentity embedding, forcing K=Q. Parallel
x=1 is excluded by the actual positive-null compactness argument. Thus
every finite, parallel and dotted x branch is covered.

## Simultaneous payments and cost

Two selected half-weight corners already pay1. With exactly one selected
corner, any unchanged positive corner pays at least1; the all-other-zero
branch was excluded above. The earlier A20/D20/B20 and balanced-quarter
changes have complete white graphs different from B16+2A1. They cannot
silently reduce a second positive angle in this quad to a fraction below1.

The inherited high-label triangle theorem supplies an unchanged entire
I2(m>=6) corner of weight2, which cannot be a whole B19+A1 vertex. The
nonnegative larger-polygon bound is untouched. This gives the required
consumer retention, including simultaneous selected angles.

Finally the B19 path has35 distance-at-most-two pairs, and its A1 has
none. The existing balanced-quarter recipient is a different marked
angle and adds1/4. Subtracting the new1/2 gives139/4. This is a verified
local whole-forest cost, not a statement that all remaining forests fit
35. Whole B20 and other unresolved D/E/high-dihedral families remain
explicit limitations. No finite-volume or constructive H9 claim is
associated with this arithmetic charge theorem.

The inherited compactness/simplicity, acute projection, finite Coxeter
classification, arithmetic cyclic-field/embedding criterion and high-label
triangle result remain proof inputs. This review independently checked
the new deductions from those inputs, without rerunning any historical
classification or claiming a new field list.
