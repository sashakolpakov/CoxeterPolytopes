# Vinberg payment for quadrilaterals with four finite corners

**Claim.** Let \(P\subset\mathbb H^n\) be a finite-sided, finite-volume Coxeter polytope which is simple at every finite and ideal vertex, with \(n\ge18\). If a quadrilateral two-face \(F\) has four finite vertices, then at least one of its corners has positive weight under Vinberg's original finite-vertex rule. In particular its total corner payment is at least one.

This is a transfer of the published diagram argument, not a new classification of its finite cases. The audit below retains Vinberg's finite Coxeter diagram classifications and local determinant calculations as published inputs and checks precisely the hypotheses that used compactness. It does not apply the compact theorem directly to a noncompact polytope.

The primary source is Vinberg, *Trudy Moskovskogo Matematicheskogo Obshchestva* **47** (1984), 68–102, especially §3 Proposition 8, §6.3, and §6.5–11, printed pp.82–91. [Original paper](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt), [bibliographic record](https://www.mathnet.ru/eng/mmo444).

The [simple-cusp affine exclusion](fresh_simple_cusp_affine_exclusion_20261017.md) is the other geometric input: every connected affine subdiagram of \(P\) has exactly \(n\) nodes. In particular every subdiagram on fewer than \(n\) nodes is affine-free. An order-\(n\) affine diagram omitting a supporting white wall is **not** excluded by this fact. No step below assumes otherwise.

## Geometric properties of the actual compact face star

Write \(W\) for the \(n-2\) supporting white walls of \(F\), and write the four black walls as opposite pairs \(\{u_1,u_2\}\), \(\{v_1,v_2\}\). The full star has \(n+2\) roots. Every \(W+u_i+v_j\) is positive definite, since it is the wall diagram at an actual finite vertex. Also \(W\) is positive definite. These assertions do not require compactness of the ambient \(P\).

Both \(W+u_1+u_2\) and \(W+v_1+v_2\) are indefinite. Here is a direct argument that deals with the potential ideal case. If the first Gram matrix were positive definite, its inverse would be entrywise nonnegative. Project a strict timelike interior point \(h\) onto the common orthogonal complement of its normals. The projection adds a nonnegative combination of those normals; its norm becomes more negative and its products with all other walls remain strictly negative. This produces an actual point of \(F\) on both opposite sides, which is impossible.

If instead that Gram matrix were positive semidefinite and singular, eliminate \(W\). Each of the two residual diagonals is positive, and the residual off-diagonal is nonpositive. Singularity gives a kernel vector with both black coordinates strictly positive. Lifting it through the entrywise nonnegative inverse of \(W\) gives nonnegative white coordinates. The resulting combination \(z\) of actual normals is null and nonzero: \((z,h)<0\) for a strict center \(h\). It is orthogonal to \(W\) and both opposite roots, and has nonpositive product with every other wall. Thus it gives an ideal point of \(F\), contradicting compactness of this four-finite-vertex face. Negative index at most one leaves the asserted indefinite alternative.

Every indefinite subdiagram of the star contains an opposite black pair, since otherwise it lies in one of the four positive-definite corner diagrams. Hence all four black roots belong to the unique indefinite connected component \(H\).

The projected side normals of a nondegenerate compact quadrilateral span the three-dimensional Lorentz space \(W^\perp\). Consequently the full star has rank \(n+1\) and corank one. The component \(H\) has corank one as well; its other components are positive definite by the orthogonality lemma below. In particular \(\det H=0\), just as required in §6.10.

The unique dependence has nonzero black coefficients, with the two opposite \(u\)-coefficients of one sign and the two opposite \(v\)-coefficients of the other. One way to see this is to normalize the projected side normals against an interior point of \(F\). Their spatial parts are the vertices of the polar convex quadrilateral. Its unique affine dependence has alternating signs around the quadrilateral. This also proves that no black coefficient vanishes.

## Two replacements for blanket affine exclusion

First, if an actual subdiagram \(A\) is orthogonal to an indefinite subdiagram \(B\), then \(A\) is positive definite. The span of \(B\) contains a timelike vector, whose orthogonal complement is positive definite. Thus the Gram matrix of \(A\) is positive semidefinite. If it were singular, a connected singular Coxeter component would have a positive kernel vector. The associated combination of its normals cannot vanish because it pairs strictly negatively with an interior point. This is impossible in the positive-definite orthogonal complement. This argument works for arbitrarily large \(A\), and needs no affine-order cutoff.

Second, an indefinite subdiagram with at most \(n\) nodes is affine-free in the present simple finite-volume setting. A proper connected affine subdiagram would have fewer than \(n\) nodes; an affine subdiagram of order \(n\) would be the whole diagram and could not be indefinite. The same assertion applies componentwise when appropriate.

These two replacements are stronger than merely observing negative index one: an abstract block-diagonal matrix could have both an indefinite and a singular positive-semidefinite block. The strict center of the actual chamber rules out the singular block here.

## Preliminary diagram decomposition

Use the decomposition \(H=U\cup V\cup W_{\rm shared}\) from §6.5: \(U\) contains the two \(u\)-roots and the white components unable to reach a \(v\)-root by a white path; \(V\) is defined symmetrically. The notation \(W_{\rm shared}\) here distinguishes Vinberg's shared partial diagram from the whole supporting diagram \(W\) above.

Assume for contradiction that all four corners have weight zero. The shortest white path between the two black pairs has length at least two, so \(U,V\) are orthogonal. There is at least one shared white node. Since each of \(U,V\) contains two black nodes and \(|H|\le n+2\),

\[
|U|,|V|\le n-1,\qquad |H-U|,|H-V|\le n.
\]

The finite corners and the opposite-pair indefiniteness proved above give the connected indefinite diagrams \(H-U,H-V\), exactly as in §6.5 Lemma 4. They are therefore affine-free by the second replacement. The diagrams \(U,V\) are affine-free by their strict size bounds. Thus §6.6 Lemma 5 still gives an elliptic member of the orthogonal pair \(U,V\), and Lemma 12 still supplies the required cyclic Lannér core whenever the corresponding elliptic member has no high bond.

All the auxiliary \(T_{ij}\) are actual finite corner components. Therefore Lemmas 1–3 and 6–11, including the end-segment conditions, the four exceptional \(E_8\) marked positions, the tree property of the shared white paths, and the classical-versus-\(E_8\) alternatives, apply to these \(T_{ij}\) literally. None of these finite-corner statements is asserted for an affine corner.

Lemma 7 is needed before bounding the order of \(U\) or \(V\): if either contains a high bond, that whole diagram is connected. Thus an **elliptic** \(U\) or \(V\) with a high bond is a connected finite diagram, not an arbitrarily large forest containing a small exceptional component. Its high bond necessarily has a finite Coxeter label at least five (multiplicity at least three in Vinberg's convention), and the finite classification gives order at most four.

Among connected affine diagrams, a **finite** bond of label at least five occurs only in the three-node affine \(G_2\) diagram, with label six. Hence a connected diagram of order at least five containing such a finite high bond is indefinite, at any order. This statement does not concern the infinite bond of affine \(A_1\), and does not use a global affine prohibition.

## Audit of the five remaining cases

The following table records the affine-sensitive parts of §6.7–11. The published local finite-diagram and determinant conclusions are retained. The cardinalities are deliberately generous carrier bounds, rather than counts of the minimal offending affine diagram.

| Part of the published proof | Justification available here |
|---|---|
| §6.7: \(U\) elliptic with a high bond, \(V\) indefinite | \(|U|\le4\). Connected high-bond extensions either become indefinite within a small extension or belong to the explicitly displayed exceptional configurations. The latter fit inside \(U\), one \(E_8\) corner component and the remaining black node, at most 13 nodes. The displayed affine obstructions have at most nine nodes. Pieces orthogonal to \(V\) are elliptic by the first replacement. |
| §6.8: \(U\) elliptic without a high bond, \(V\) indefinite | The cyclic Lannér core has at most five nodes. Its local extensions, including the two three-node end pieces, fit in a carrier of at most 13 nodes. The explicit affine \(B_4/F_4\) exclusions are small. Any larger piece shown orthogonal to \(V\) is already elliptic by the first replacement. |
| §6.9: both \(U,V\) elliptic with high bonds | \(|U|,|V|\le4\). The exceptional \(E_8\)-corner configurations fit inside \(U\cup V\cup T_{ij}\), at most 16 nodes. The displayed affine obstructions have at most nine nodes. For the long classical case, see the separate argument immediately below. |
| §6.10: only one of the two elliptic diagrams has a high bond | The exceptional configurations are contained in \(U\cup T_{ij}\cup L\), at most \(4+8+5=17\) nodes. Other contradictions use small high-bond extensions orthogonal to the cyclic Lannér core, or the negative-index and determinant tests. The quadratic determinant equation is nontrivial by the argument below. |
| §6.11: neither elliptic diagram has a high bond | The two cyclic Lannér diagrams and all their overlap/connecting edges use at most ten nodes. The implicit affine exclusions used to determine their permitted overlap or connection are within this carrier. The later larger end-piece contradictions use the already known ellipticity of \(U,V\) and the finite \(T_{ij}\) classification; they do not apply a blanket affine exclusion to those larger diagrams. |

Every finite carrier in this table has fewer than \(n\) nodes when \(n\ge18\). Its affine exclusions therefore remain valid.

For precision in \S6.11, the small-carrier assertion concerns restrictions on the induced union of the two cyclic cores and their connecting edges. It is not a separate assertion that arbitrary white branches have been excluded by an affine test on that union. The subsequent determination of the shared paths uses the already positive-definite \(T_{ij}\) and Lemma 8: a cross-distance-two shared diagram is linear, and at the other short cross-distances a nonlinear shared diagram requires the finite \(E_8\) alternative. The published end-piece arguments then use Lemmas 10--11 and ellipticity of \(U,V\). Thus extending a displayed core pattern by additional white nodes does not require a new blanket affine exclusion; the pertinent larger diagrams are already finite corner diagrams. This distinguishes the core-overlap affine tests from the later finite classification of complete shared paths.

There is one additional size point in the long classical branch of §6.9. Its (O\(U\)) and (O\(V\)) attach every white neighbor of the two black roots to \(U,V\). Each black root has at most three white neighbors, because all these neighbors occur in an actual positive-definite corner diagram. Thus

\[
|O\(U\)|,|O\(V\)|\le4+2\cdot3=10<n.
\]

They are orthogonal when the shortest white cross-path has length at least four, so one is elliptic. After the long corner component is split into its two end pieces, the high-bond pieces of order at least five are indefinite by the finite/affine classification stated above. Their orthogonal complements are elliptic by the first replacement, regardless of their order. This supplies the ellipticity needed for the final \(H_4\)-piece and negative-index contradiction in that branch. No order-\(n\) affine end piece has been silently discarded.

## Nonzero determinant polynomial in §6.10

Let \(e\) be the high bond of \(U\). All the other normalized Gram entries belong to \(\mathbb Q(\sqrt2)\): the elliptic \(U\) has only this one high bond, the elliptic \(V\) has none, and every other relevant bond belongs to the finite corner diagrams covered by Lemma 6. Varying the entry on \(e\) gives a determinant polynomial of degree at most two. In the original diagonal-one normalization its quadratic coefficient is minus the determinant obtained by deleting the endpoints of \(e\). Its positivity deserves an explicit check in the present setting.

In the unique dependence of \(H\), every vertex of \(U\) has a nonzero coefficient of the same sign as its two black \(u\)-coefficients. Indeed, a private white component has contacts only with the \(u\)-pair and its own white nodes. Solving its dependence equation using its entrywise nonnegative positive-definite inverse gives that common sign, strictly at every node because the component is connected and attached to the black pair.

Deleting an endpoint of \(e\) therefore destroys the only linear dependence. The normals remaining after both endpoints are deleted are independent, and their span still contains the indefinite diagram \(H-U\). A subspace containing a timelike vector has nondegenerate Lorentz restriction. Its Gram determinant is consequently strictly negative. Thus the quadratic coefficient is strictly positive (four times as large if the variable is \(\cos(\pi/m)\) in diagonal-two normalization), including when the high bond has private white endpoints.

The published algebraic angle argument and the subsequent negative-index contradictions can therefore be used unchanged. This check relies on the actual compact quadrilateral and its full star; it would not be justified by deleting an arbitrary white root and pretending the smaller diagram were another face star.

## Conclusion and review scope

Under the assumption of four zero-payment finite corners, the preliminary decomposition and all five published cases retain their required hypotheses. Each case produces the same contradiction. This proves the claim for four-finite-vertex quadrilaterals in the stated all-simple finite-volume class, \(n\ge18\).

This note is an independent mathematical/source audit, not a numerical graph census or an independent rederivation of every determinant table in the original paper. Separate source and mathematical reads by the root and frontier agents pass the geometric replacements, Lemma 7 connectedness, carrier-size controls, §6.11 locality, and the nonzero quadratic coefficient. No remaining compactness gap was identified in the stated transfer. These internal reviews do not claim to duplicate the published local tables or constitute external refereeing or formal verification. This proof note is final and frozen.

The other mixed-corner quadrilateral branches require their separate proofs. This note alone gives no new ambient dimension bound, and changes none of the frozen 61/127/163 proofs or arithmetic certificates.
