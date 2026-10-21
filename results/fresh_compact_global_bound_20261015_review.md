# Independent review of the compact cubic-seven dimension-20 candidate

**PASS within the stated compact arithmetic scope.** The [standalone checker](../scripts/fresh_compact_global_bound_20261015_review.py) and [certificate](fresh_compact_global_bound_20261015_review.json) independently reconstruct the finite graph reduction, polynomial rank identities, exceptional corner weights and dimension inequality in the [primary proof](fresh_compact_global_bound_20261015.md). Generation and saved replay pass. This review finds no remaining mathematical gap in that proof, subject to its explicit published inputs.

The conclusion concerns compact hyperbolic Coxeter polytopes whose arithmetic ground field is \(\mathbb Q(\cos(2\pi/7))\): their dimension is at most 20. It does not improve the unconditional compact bound, remove either quadratic arithmetic field, or claim priority. The previous [dimension-28 fallback](fresh_arithmetic_highdim_20261015.md) and its [independent review](fresh_arithmetic_highdim_20261015_review.md) remain frozen.

The new checker imports no producer code and uses no SymPy. A separate Fraction polynomial implementation expands both literal Gram determinants by sparse Laplace recursion, reduces powers using the cubic relation, and verifies every coefficient identity needed by the proof. It enumerates complete binary graphs instead of the producer's connected-shape/cross-edge loops. The earlier field classification, standard compact Coxeter face facts and the published arithmetic criterion are declared dependencies rather than claimed as new independent classification results.

## Domain and geometric completeness

The global compactness assumption is essential and is used correctly. It excludes affine subdiagrams and supplies the actual-face correspondence for positive definite wall subsets. At a quadrilateral two-face, each adjacent black pair with the complete white support is spherical, whereas each opposite extension is indefinite. Its minimal indefinite diagram must contain both opposite black nodes, since white support plus either single black node occurs at a finite vertex. The inherited order bound therefore leaves exactly orders two and three.

The distinguished high bond in an order-three block does not require trusting the auxiliary seven-triangle enumeration. The inherited finite alphabet is \(2,3,4,6,7,14\). If a Lannér triangle had no 7 or 14 label, its determinant would lie in \(\mathbb Q(\sqrt2,\sqrt3)\); it could not generate a cubic ground field. This contradicts the same determinant-field proposition used in the fallback. Thus every such triangle has a 7/14 bond. A bond of label at least six is isolated in any spherical diagram containing it.

The review enumerates all 35,904 binary graphs on the five possible minimal-block vertex configurations, including shared white nodes and every cross-bond choice. For each graph it checks connectedness of both minimal blocks, finite-corner forest constraints, isolation of each possible distinguished high edge, and absence of a marked path of length at most two. Exactly the producer's nine supports survive. They all have two distinct white nodes joined by an edge, with both high edges on opposite black pairs. In particular neither a dotted minimal pair nor a high black-white bond has been omitted.

The extension from these initial nodes to **all** white contacts is valid. Project to the orthogonal complement of the two positive definite opposite black blocks. A white node contacting either block has residual square below \(-6/19\). A 14 block or a dense contact to both roots of a 7 block makes one square below \(-4\). The cross product between white nodes attached to different blocks remains their old product, because the absence of short corners forbids cross contacts. Its square is less than one. The resulting two-negative-plane determinant is then greater than \(5/19\), a contradiction. This excludes 14 and dense attachments at every white node, not only the initially chosen ones.

Let \(U,V\) collect all remaining white neighbors of the respective opposite black pairs. Every \(U\)-\(V\) edge is forced nonright: otherwise the two negative residual vectors would be orthogonal. The white diagram is a forest, so one part is a singleton. A black endpoint cannot have two neighbors in the other part, since this would make a cycle in its finite corner diagram. Consequently the complete touched white set has sizes \(1+1\) or \(1+2\). This rules out hidden larger support patterns.

Every off-core white component attaches at most once to the connected white core; two attachments would produce a white cycle. Its connected edges have labels three or four, since label five is unavailable and a larger label cannot be attached to another edge in a finite diagram. Tree elimination therefore has rational diagonal corrections. A nonzero branch correction is strictly positive, and all reduced white diagonals remain positive and at most one.

## Rank, rationality and the entire-component condition

The complete star has \(|W|+4\) roots in a space of dimension at most \(|W|+3\): the orthogonal complement of the white support is Lorentz three-space. After elimination of off-core white roots, the six- or seven-node matrix must therefore be singular. Only this upper rank bound is needed. A smaller local diagram with Lorentz signature but excessive rank is not a realizable two-face star.

The independently expanded six-node determinant is a nonzero rational quadratic in \(c^2\), with leading coefficient the positive reduced white determinant. Since \(c^2\) has degree three, this cannot vanish. This correctly removes the apparent path \(7\!-!3\!-!3\!-!3\!-!7\), which passes some local tests but fails the complete-star rank condition.

In the seven-node case the rationality claim also survives scrutiny. Squared bonds and reduced diagonals are rational. The one potentially irrational coefficient is \(L=\delta\epsilon\beta\gamma\in\mathbb Q(\sqrt2)\). The actual five-cycle has product a nonzero rational multiple of \(Lc\). The arithmetic cyclic-product condition and \(c\in K\) put \(L\) in \(K\), whose intersection with the quadratic field is rational.

The separate Fraction determinant expansion gives the linear coefficient of the cubic remainder as \((\det H-8L)/8\). Hence singularity forces \(\det H=8L\). Positive definiteness makes the determinant strictly increasing in each white diagonal, and therefore

\[
\det H\le1-\delta^2-\epsilon^2\le\tfrac12\le8L.
\]

Equality requires all three diagonals to be one and \(\delta=\epsilon=\beta=\gamma=1/2\). No nonzero eliminated branch can remain. The other coefficient then gives \(\alpha=1/2\). These claims are checked coefficientwise, not by numerical substitution at the cubic root.

The remaining four corners have two types \(D_5\) and two types \(A_4+A_1\). Independent rational LDL tests certify all four. The marked \(D_5\) pairs are precisely the two far-long-leaf/short-leaf pairs at distance three. All black-white contacts have already been exhausted by \(U,V\), and all external white branches were removed by equality. Thus each weighted \(D_5\) is an **entire connected component** of its full finite vertex diagram. This is the necessary justification for the half-weights; an arbitrary embedded \(D_5\) would not suffice.

## Global weight accounting and arithmetic scope

The triangle and quadrilateral lower weights are correct: a triangle contains a high rank-two corner receiving two units; a quadrilateral either has a short corner receiving one or the classified pair of \(D_5\) half-weights. Larger polygons need only nonnegative weights.

The review reconstructs 122 finite diagrams. The all-rank proof follows directly from edge and neighbor-pair counts, not from that finite range. Every connected component has deficit at least two from \(2r\), except \(D_5\), whose deficit is one. At total rank at least six, a single component cannot be the exceptional \(D_5\); two or more components have total deficit at least two. Hence the vertex weight is at most \(2n-2\).

The strict h-vector lower bound is applied with the correct parity. At dimension 21 it is strictly 40, matching the upper bound 40 and giving a contradiction. The independently reproduced gaps for all larger odd/even dimensions are \(u^2/2+5u\) and \(u^2/2+11u/2+2\). Strictness handles the zero constant in the odd case. The resulting bound is 20.

The arithmetic criterion and weighted-face/h-vector framework are the stated imports from [Vinberg's original paper, §§1–3 and 7](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt). This review is a reproducible mathematical and exact-arithmetic audit, not a formal proof-assistant verification or an exhaustive search for prior equivalent results.

```sh
python3 scripts/fresh_compact_global_bound_20261015_review.py --verify
```
