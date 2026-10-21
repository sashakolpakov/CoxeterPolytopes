# A distinguished-bond obstruction in eight arithmetic fields

**Research proof candidate.** A compact arithmetic hyperbolic Coxeter polytope has dimension at most **20** if its ground field K satisfies

\[
 K\not\cong\text{a subfield of }F=\mathbb Q(\sqrt2,\sqrt5),
 \qquad \mathbb Q(\cos(2\pi/7))\not\subseteq K. \tag{1}
\]

Among Vinberg's twelve possible compact fields in dimensions at least 14, this covers eight fields:

\[
\mathbb Q(\sqrt3),\quad\mathbb Q(\sqrt6),\quad
\mathbb Q(\sqrt2,\sqrt3),\quad K_9,K_{11},K_{15},K_{16},K_{20},
\qquad K_m=\mathbb Q(\cos(2\pi/m)).
\]

The proof below is independent of the large Q(√5) attachment queues and does not assert a dimension-20 bound for all arithmetic chambers. The [standalone checker](../scripts/fresh_arithmetic_high_label_20261016.py) and [certificate](fresh_arithmetic_high_label_20261016.json) verify the finite support reduction, scalar identities, field-group facts, graph counts and strict dimension threshold. Generation and fresh saved replay pass. Independent review and external validation are separate obligations; no novelty or current-best claim is made.

## Published inputs and exact field scope

We use whole-chamber compactness, the positive-definite face correspondence for acute chambers, finite Coxeter diagrams, and the symmetric nonnegative h-vector of a simple polytope. Compact Coxeter chambers are simple and have no affine principal diagrams. A minimal indefinite principal diagram is Lannér. Its normalized determinant generates the arithmetic ground field. These are the relevant inputs from [Vinberg, §§1–3 and Proposition 17](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).

For comparison, [Nikulin, introduction and §§2–3.2](https://arxiv.org/pdf/0708.3991) states the dimension-22 field restriction and the larger dimension-14 list, explains cyclic products, and reproduces the arithmetic simplex/triangle field data. This proof does **not** assume that a field restriction valid from dimension 22 remains valid in dimension 21, and does not need to import the full triangle catalogue.

Normalize all normals to square one. A finite bond m then has product \(-\cos(\pi/m)\). The length-two cyclic product gives \(4\cos^2(\pi/m)\in K\). A normalized determinant is a sum of products of cyclic products and rational constants, and hence is in K. For a Lannér determinant, its negative physical value and positive values under all other ground-field embeddings force it to generate K; field **equality**, not merely containment, is essential here.

## 1. What (1) forces on every Lannér diagram

**No Lannér diagram has order at least four.** Such a diagram is connected and all its proper principal subsets are positive definite. A bond m≥6 has a neighboring nonright bond somewhere along the connected diagram. The resulting proper three-node diagram has determinant at most

\[
1-\cos^2(\pi/6)-\cos^2(\pi/3)=0;
\]

any third bond only decreases the determinant. Dotted or parallel bonds already violate a proper two-node test. Thus every bond has label at most 5. All its Gram entries, and its determinant, lie in F. Since that determinant generates K, this contradicts (1). This reasoning does not infer arithmeticity from a low-degree determinant; it uses the arithmeticity of the hypothetical ambient chamber.

Likewise, a Lannér triangle with all labels at most 5 has determinant in F and is impossible. A label 7 is impossible anywhere, since its length-two cyclic product generates K₇. Consequently, every Lannér triangle satisfies one of the following:

1. it contains a bond of label at least 9; or
2. all its labels are at most 8, it contains a bond 6 or 8, and **all three bonds are nonright**.

For the last assertion, a right triangle with labels from \(\{2,3,4,5,6,8\}\) has determinant \(1-\cos^2(\pi/a)-\cos^2(\pi/b)\in F\), again impossible. In particular, a label 6 or 8 in the second case has a nonright contact at each endpoint.

We call either chosen bond a distinguished bond. Every distinguished bond has label at least 6, so it is an **entire rank-two component in every finite diagram containing it**: another nonright incident bond would violate the same three-node positive-definiteness inequality.

## 2. The two relevant Schur costs are strictly larger than two

Let an opposite finite pair have Gram

\[
B_c=\begin{pmatrix}1&-c\\-c&1\end{pmatrix},\quad
B_c^{-1}=\frac1{1-c^2}\begin{pmatrix}1&c\\c&1\end{pmatrix},
\quad c=\cos(\pi/m).
\]

A third normal with nonpositive pairings \((-p,-q)\) has projection cost

\[
\frac{p^2+q^2+2cpq}{1-c^2}. \tag{2}
\]

Every nonright contact has magnitude at least 1/2.

If m≥9, at least one contact suffices. Put \(z=2\cos(2\pi/9)>1\). It satisfies \(z^3-3z+1=0\); that polynomial is strictly increasing for z>1 and equals −1/8 at 3/2. Hence z>3/2 and \(\cos^2(\pi/9)>7/8\). Formula (2) is therefore strictly larger than

\[
\frac{1}{4(1-7/8)}=2. \tag{3}
\]

If m≥6 and both contacts are nonright, (2) is at least

\[
\frac1{2(1-c)}>2, \tag{4}
\]

because \(c\ge\cos(\pi/6)=\sqrt3/2>3/4\). Thus every distinguished bond costs more than two, including the dense label-6/8 alternative. The latter is necessary: the arithmetic triangle (3,3,8) has field K₁₆ and would defeat a restriction to labels at least 9 alone.

## 3. Every quadrilateral has a marked path of length at most two

Let W be the white support of a compact quadrilateral two-face and let a,b,c,d be its black side nodes in cyclic order. The four W-plus-adjacent-pair diagrams are positive definite. Both opposite extensions are indefinite: positive definiteness would give an actual intersection of opposite sides, while a singular positive-semidefinite subset would give an affine diagram in the compact chamber.

Choose minimal indefinite diagrams L⊆W+{a,c} and M⊆W+{b,d}. Each contains its two opposite black nodes, and by Section 1 has order two or three. Suppose that no adjacent black pair has distance at most two in its finite corner diagram.

All adjacent black contacts are right. If L or M has order two, every cross contact between the two selected diagrams is forbidden: a black-to-black contact gives a length-one marked path, and a black-to-white contact gives a length-two path using a nonright bond of the connected triangle. The two indefinite diagrams would then be orthogonal, giving two negative directions. Therefore both are triangles.

Their white nodes u,v are distinct, since a shared white node would give a length-two path between some adjacent black pair. The same path argument forbids all cross black-to-white and black-to-black contacts. The only possible connection between L and M is u–v, and it must be nonright to avoid orthogonal indefinite blocks.

A distinguished bond incident to u or v would share its white endpoint with u–v inside a finite corner. This contradicts rank-two isolation. Thus the distinguished bonds are exactly the opposite black pairs a–c and b–d. Their Gram is a positive-definite direct sum of two rank-two blocks.

Eliminate that four-node block in the principal six-node Gram of a,b,c,d,u,v. By (3) or (4), both residual white diagonals are strictly smaller than −1. The off-diagonal is unchanged: u and v have no contacts to the other one's opposite pair. Its square is strictly smaller than one, because the white diagram is positive definite. Consequently the residual 2×2 matrix has negative diagonal and strictly positive determinant. It is negative definite, contradicting the ambient Lorentz negative index one.

This excludes the supposed quadrilateral. Arbitrary other white nodes cannot fix a principal submatrix with two negative directions. **No complete-star rank assumption, finite white-rank cutoff, sealing assumption, or spectator enumeration is needed for this local lemma.**

The checker independently lists all binary supports of order-two/three minimal pairs, their shared-white possibilities, every cross-contact subset and a chosen distinguished bond on each triangle. Relaxing finite corners to forests with distinguished-bond isolation leaves exactly nine supports: the three possible nonempty white attachment patterns on each of the two opposite pairs. The calculation verifies the incidence argument; it is not a substitute for the continuous cost inequalities.

## 4. Weights and the all-rank contradiction

For a planar angle, mark the two incident facet nodes defining its two-face in the complete finite vertex diagram. Give it:

- one unit when the marked nodes have graph distance at most two;
- one further unit when they form an entire I₂(m) component with m≥6;
- zero otherwise.

The second term is only on rank-two components and the weights are nonnegative. Every triangular two-face contains a Lannér diagram containing all three black nodes; it is exactly that triangle. Its distinguished bond is an isolated I₂ component at the corresponding corner, giving total payment at least two. Section 3 gives at least one unit on every quadrilateral. Larger polygons have nonnegative payment, at least their required value 5−k.

For a connected finite diagram of rank r, the number of pairs at distance at most two equals

\[
(r-1)+\sum_v\binom{\deg v}{2}.
\]

It is 2r−3 on a nontrivial path, 2r−2 on a D/E tree, and zero on A₁. The exceptional finite F/H diagrams are paths. The added dihedral unit makes the rank-two total two, still 2r−2. Therefore every connected component has total at most 2r−2, and every nonempty rank-n vertex forest has total at most 2n−2. The all-rank statement follows from these degree formulas; the 126 explicit graph records are checks of the formulas, not a rank extrapolation.

For a simple n-polytope the total face weight is at least \(5f_2-\binom n2 f_0\). Symmetry and nonnegativity of its h-vector give the strict average lower bounds

\[
\frac{n(n-6)}8\quad(n\text{ even}),\qquad
\frac{(n-1)(n-5)}8\quad(n\text{ odd}).
\]

Strictness follows from h₀=hₙ=1 and the larger extreme paired coefficient. At n=21 the lower bound is strictly greater than 40 while the vertex upper bound is 40. For n=21+2u and n=22+2u, the non-strict lower-minus-upper differences are respectively \(u^2/2+5u\) and \(u^2/2+11u/2+2\). Thus all n≥21 are excluded.

## 5. Eight fields and the precise remaining dimension-21 branches

The only subfields of F are Q, Q(√2), Q(√5), Q(√10), and F itself. The fields Q(√3), Q(√6), and Q(√2,√3) are not among them. K₉ and K₁₁ have degrees three and five. K₁₅,K₁₆,K₂₀ are cyclic quartic fields, whereas F is biquadratic. These establish the first condition of (1).

All eight except K₉ have degree not divisible by three, so cannot contain K₇. For K₉, the cyclotomic fields of conductors 7 and 9 have product degree \(\varphi(63)=\varphi(7)\varphi(9)=36\), hence intersect only in Q. This establishes the second condition. The checker reconstructs the real cyclotomic quotient groups from units modulo m and checks the degrees/orders and squareclasses explicitly.

The exact compact twelve-field list at dimension at least 14 is

\[
\mathbb Q(\sqrt2),\mathbb Q(\sqrt3),\mathbb Q(\sqrt5),\mathbb Q(\sqrt6),
\mathbb Q(\sqrt2,\sqrt3),\mathbb Q(\sqrt2,\sqrt5),
K_7,K_9,K_{11},K_{15},K_{16},K_{20}.
\]

The frozen independently reviewed [Q(√2)](fresh_sqrt2_quad_exhaustion_20261015_review.md) and [K₇](fresh_compact_global_bound_20261015_review.md) arguments already exclude dimension 21 in those fields. If the new proof passes independent review, the only remaining compact dimension-21 fields are **Q(√5)** and **Q(√2,√5)**. The current complete [Q(√5) bound](fresh_sqrt5_h3_hub_coverage_20261016.md) is only 21 and does not exclude dimension 21. In F, the arithmetic triangle (2,5,8), present in the cited triangle list, lies outside the distinguished-bond mechanism: the 8 bond has only one nonright attachment, and the squared 5 and 8 contacts already generate F. Neither unresolved branch is silently removed by the dimension-22 three-field restriction.

This is an explicit new proof candidate and branch reduction, not a claim that compact dimension 21 has been excluded or that no prior source proves the same subclass result. A bounded primary-literature check did not establish priority.

```sh
python3 scripts/fresh_arithmetic_high_label_20261016.py
python3 scripts/fresh_arithmetic_high_label_20261016.py --verify
```
