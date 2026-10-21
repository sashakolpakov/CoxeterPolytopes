# All-real obstruction to the one-contact A4-to-H4 section family

**Result.** No hyperplane section of the saved Borcherds H21 chamber can
be a finite-volume Coxeter H20 polytope if it has an actual H4 face obtained
from four norm-2 A4 walls by changing only one endpoint pairing. This
excludes every real position and both branches in all **3,320** such
families attached to the **1,660** A4 components in the saved source-edge
census. Section 6 extends the obstruction to all **4,230** one-contact
A3-to-H3 families attached to its **2,115** A3 components. It supersedes
the earlier nine exact failures and 108 numerical
samples for this family. It does not exclude other local section diagrams,
multiple-contact modifications, or other ambient polytopes.

The argument uses compactness and rational arithmetic, not a finite
parameter sample. The [exact checker](../scripts/fresh_h20_h4_face_closure_20261014.py)
and [certificate](fresh_h20_h4_face_closure_20261014.json) verify its finite
source and Gram-matrix premises. They do not mechanically formalize the
standard geometric and arithmetic theorems used below.
The [independent SymPy and proof review](fresh_h20_h4_face_closure_20261014_review.md)
passes on the frozen inputs, including the H3 extension.

## 1. Source data and the section

Use the integral Lorentzian lattice L=Z^22 with form Q of signature (21,1)
and the 210 outward roots r_i in the
[source packet](../data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json).
The exact checker verifies that every root reflection preserves L:
2Qr_i/(r_i,r_i) has integral coordinates. There are two root norms, 2 and 4,
with precisely these off-diagonal pairing pools:

| Root norms | Possible pairings |
|---|---|
| 2,2 | −1, 0 |
| 2,4 | −2, 0 |
| 4,4 | −8, −6, −4, −2 |

Let C be a norm-2 A4 chain and a one of its endpoints. Scale the section
normal e so that (e,r_a)=1 and (e,r_i)=0 for the other three nodes of C.
For its projected chain to be H4, the squared normal norm must be

\[
q=(e,e)=\frac{1+\sqrt5}{4}.
\]

Indeed the changed endpoint norm is 2−1/q, and the adjacent squared
cosine is 1/[2(2−1/q)]=(3+√5)/8=cos²(π/5). The other chain labels stay 3.
Write s_i=r_i−(e,r_i)e/q for projected normals.

In the source-edge construction, let w_a be the endpoint dual vector in
the span of the twenty edge walls. Then w_a²=4/5, and every such section
has e=w_a+z, where z is in the Lorentzian two-plane orthogonal to the
source edge walls and

\[
z^2=q-\frac45=\frac{5\sqrt5-11}{20}>0.
\]

An interior transverse edge crossing gives an actual finite H4-containing
vertex. The argument below requires only that the projected H4 is an
actual codimension-four face, not that the other sixteen edge pairings
vanish. In particular it applies to all these edge families.

## 2. The H4 face is a compact Coxeter H16 polytope

Assume for contradiction that the section P is a finite-volume Coxeter
polytope. Let F be its face cut out by the four H4 walls.

H4 has no proper connected spherical extension. At a finite point of a
facet of F, the extra wall therefore meets all four H4 walls orthogonally.
Moreover H4 is contained in no irreducible affine Coxeter diagram. The
bounded Euclidean link of an ideal vertex has only affine components, so
F has no ideal vertex. A face of a finite-volume hyperbolic polytope is
the convex hull of its finite and ideal vertices. Thus F is compact.

All its boundary walls are orthogonal to the H4 walls, so their mutual
angles inside F equal their angles in P. Hence F itself is a compact
Coxeter polytope of dimension 16. This also follows from the usual
Coxeter-face criterion; the elementary orthogonality argument specifies
the normals needed here. The spherical/affine classification and the
Coxeter-face construction are recalled in
[Felikson's research overview](https://www.maths.dur.ac.uk/users/anna.felikson/Polytopes/polytopes.html)
and [Borcherds' paper](https://arxiv.org/abs/math/9806141).

For each source wall l inducing a facet of F, set k_l=(r_l,r_a).
Orthogonality to s_a and to the other three H4 normals implies

\[
(r_l,r_i)=0\quad(i\in C\setminus\{a\}),\qquad
(e,r_l)=qk_l.
\]

Its normal inside F is therefore

\[
s_l=r_l-k_le,\quad s_l^2=r_l^2-qk_l^2,\quad
(s_l,s_m)=(r_l,r_m)-qk_lk_m.
\]

If k_l=0, this is literally the original integral root r_l. Call these
facets unmodified. By the exact source pairing pools, a modified facet
has either (r_l²,k_l)=(2,−1) or (4,−2).

## 3. Every neighbor of a modified facet is unmodified and orthogonal

For all root-type pairs involving a modified facet, the exact checker
enumerates the complete 18-row projected Gram table. Twelve rows cannot
meet at a finite point. Three rows are orthogonal pairs consisting of one
modified facet and one unmodified original root. The remaining three
have squared dihedral cosines

\[
\frac{7+\sqrt5}{22},\qquad
\frac{7+\sqrt5}{11},\qquad
\frac{29+12\sqrt5}{121}.
\]

None is a squared Coxeter cosine. Here there is no angle-label cutoff:
if cos²(π/m) belongs to Q(√5), then cos(2π/m) belongs to that quadratic
field. For m≥3 its degree is φ(m)/2, so φ(m)≤4. The resulting possibilities
are m=3,4,5,6,8,10,12; m=8 and 12 generate Q(√2) and Q(√3), respectively,
and are excluded. Including m=2, the complete list is

\[
0,\quad \frac14,\quad \frac12,\quad
\frac{3+\sqrt5}{8},\quad\frac34,\quad
\frac{5+\sqrt5}{8}.
\]

The checker compares all three exceptional values exactly against this
list. Therefore every facet adjacent to a modified facet of compact F
must be unmodified and perpendicular to it.

## 4. The forced rational compact face is impossible

If F has no modified facet, every boundary normal is an original integral
root and F is a compact Coxeter H16 polytope with rational normals.

Otherwise choose a modified facet T of F. By the preceding step, every
facet bounding T comes from an unmodified root r_m perpendicular to its
normal s_l. These roots lie in the Lorentzian span of T, without any
additional projection. Angles in T are inherited from F, so T is a
compact Coxeter H15 polytope whose boundary normals are again original
integral roots.

In either case let U be the Lorentzian linear space of this compact
polytope, of dimension 17 or 16. Compactness forces its boundary normals
to span U. Since these normals are rational, U is a rational subspace and
Λ=U∩L is a full-rank integral lattice. Each boundary reflection preserves
both L and U, hence Λ. The reflection group of the compact Coxeter
polytope is a cocompact lattice subgroup of O(Λ), and is consequently of
finite index in that arithmetic lattice.

But an indefinite rational quadratic form of rank at least five is
isotropic (Meyer's theorem). By the arithmetic compactness criterion,
O(Λ) is noncocompact. A finite-index subgroup cannot be cocompact: this is
the contradiction. The rational-field compactness criterion in precisely
this hyperbolic setting is stated in
[Belolipetsky, §2, p.5](https://arxiv.org/pdf/1506.03111).

Rational coordinates alone were not used to assert arithmeticity. The
argument explicitly verifies an invariant integral lattice and uses the
finite-volume subgroup to obtain finite index.

## 5. Finite candidate reduction and numerical cross-check

Before the uniform obstruction above was derived, a separate finite
reduction replaced arbitrary sampling of z. Compact F requires a facet
outside the initial sixteen complementary edge walls. For any such wall
l the orthogonality equations give

\[
z\cdot r_l=(q-4/5)(r_l,r_a),\qquad z^2=q-4/5.
\]

This is a line intersected with a spacelike hyperbola in the edge's
Lorentzian two-plane. The equation cannot be vacuous: a wall outside the
twenty source edge walls cannot contain the entire open source edge.
There are at most two solutions per eligible l, before checking that the
section actually crosses the open edge.

The [numerical scout](../scripts/fresh_h20_h4_face_closure_20261014_scout.py)
covered all 3,320 endpoint families in the saved 6,492 literal edge
supports. It retained 152,368 candidates after within-family duplicate
removal; every candidate had a numerical forbidden angle at an actual
neighboring finite vertex. The full
[numerical packet](fresh_h20_h4_face_closure_20261014_scout.json) preserves
the recipes and witnesses. These numerical predicates are not a proof
of exhaustive root separation or exact rejection. **The all-real
exclusion rests on §§1–4, independently of this numerical packet.**

The source-edge inventory is the saved exact corner census, not a new
claim that its literal supports are inequivalent under the full source
symmetry group. The uniform proof itself does not need an edge census.

## 6. Extension to all one-contact A3-to-H3 sections

Let three norm-2 source walls a,b,c form A3, with the changed pairing at
endpoint a. Again (e,r_a)=1, (e,r_b)=(e,r_c)=0 and e²=q. Their projections
form H3 with the label5 edge at a. If this is an actual face of a
finite-volume Coxeter H20 section, its H17 face is compact: H3, like H4,
is contained in no irreducible affine diagram.

The only proper connected spherical extension of H3 is H4. Such a
boundary facet m attaches by label3 at the far endpoint c and is
orthogonal to a and b. Put k=(r_m,r_a), d=r_m² and g=(r_m,r_c). The
orthogonality to a gives (e,r_m)=qk. The new label3 gives

\[
d-qk^2=2g^2.
\]

Since d,k,g are integers and q is irrational, k=0. The source pairing
pools then force d=2 and g=−1; the norm4 case would require an integer
square equal to 2. Thus the four original source walls form A4, with
only the a-pairing changed by e. This is an actual H4 subface excluded
by §§1–4.

Consequently all boundary facets of the compact H3 face would have to
be orthogonal to H3. That face is then itself a Coxeter H17 polytope.
The same complete modified-facet table applies. If all facets are
unmodified it is rational arithmetic and compact; if a modified facet
exists, that compact H16 facet has only unmodified rational boundary
roots. The invariant-lattice and isotropy contradiction applies again.

The exact source-edge census contains 2,115 A3 components, giving 4,230
endpoint families. All real transverse section parameters in those
families are therefore excluded as well. This extension does not use
the A4 numerical scout or impose any arithmetic hypothesis on the
original hypothetical H20 section.
