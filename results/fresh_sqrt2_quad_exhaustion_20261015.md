# A complete quadrilateral reduction over Q(√2)

**Internally independently reviewed proof candidate:** a compact arithmetic hyperbolic Coxeter polytope with ground field \(\mathbb Q(\sqrt2)\) has dimension at most **20**. The finite calculations below are complete within their stated roles. The argument for arbitrarily many external white nodes is given separately; a cutoff on those nodes is not used. No unconditional compact bound or bound for other arithmetic fields is asserted. Internal verification does not establish literature priority or replace external mathematical review.

The computation consists of a [minimal-pair union enumeration](../scripts/fresh_sqrt2_quad_exhaustion_20261015.py) with [data](fresh_sqrt2_quad_exhaustion_20261015.json), an [exact arithmetic/signature filter](../scripts/fresh_sqrt2_quad_exhaustion_20261015_cores.py) with [data](fresh_sqrt2_quad_exhaustion_20261015_cores.json), a [complete one-extra-white contact calculation](../scripts/fresh_sqrt2_quad_exhaustion_20261015_extensions.py) with [data](fresh_sqrt2_quad_exhaustion_20261015_extensions.json), and a [closure/weight checker](../scripts/fresh_sqrt2_quad_exhaustion_20261015_closure.py) with [data](fresh_sqrt2_quad_exhaustion_20261015_closure.json). The [independent global review](fresh_sqrt2_quad_exhaustion_20261015_review.md) reconstructs the union enumeration, marked orbit coverage, exact core/profile arithmetic and closure identities without producer imports, and separately audits the all-size proof.

## Published inputs and the finite reduction

We use the usual compact Coxeter facts: a vertex diagram is spherical; every elliptic subset corresponds to an intersecting collection of facets; affine principal subdiagrams are absent; and a minimal indefinite subdiagram is Lannér. Arithmetic cyclic products lie in the ground field and are algebraic integers; every other ground-field embedding makes the Gram matrix positive semidefinite. These are inputs, not new conclusions of the search. [Vinberg, §§1–3 and 7](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).

The [complete Q(√2) Lannér catalogue](fresh_sqrt2_lanner_catalogue_20261015.md) proves that the finite alphabet is \(2,3,4,6,8\). Its minimal indefinite pieces are a dotted pair, nine triangle types, or a chordless C4/C5 with exactly one 4-bond and all other bonds 3. The nine unordered triangle types are

\[
(2,3,8),(2,4,8),(2,6,8),(2,8,8),(3,3,4),
(3,8,8),(4,4,4),(4,6,6),(4,8,8).
\]

The dotted distance remains an arbitrary permitted real arithmetic parameter. It is not replaced by a finite list.

Let \(W\) be the **full** white support of a quadrilateral two-face, and label its black sides cyclically \(0,1,2,3\). Every \(W+\{a,b\}\), with \(a\in\{0,2\}\), \(b\in\{1,3\}\), is spherical. The opposite extensions \(W+\{0,2\}\) and \(W+\{1,3\}\) are indefinite, and contain minimal Lannér subdiagrams \(L,M\) containing their respective two black nodes.

If any marked adjacent black pair has graph distance at most two in its full finite corner, the basic weight below already pays for this quadrilateral. Assume henceforth that no such short marked path exists. Enumerate \(L\cup M\), not all of \(W\). Each piece uses at most three white nodes, so this union has at most six whites. Both shared white nodes and all additional contacts on the union must be considered.

The enumerator includes:

- all 51 marked piece representatives, fixing the ordered two black labels and quotienting only white names;
- all 6,010 shared-white identifications, with injective identifications and canonical names for new whites;
- every unspecified union contact in the complete finite alphabet;
- all four exact finite-corner tests and absence of a marked path of length one or two.

All pairs except the two opposite black pairs lie in at least one finite corner and hence really do have finite labels. The opposite pair in each piece is already specified by that piece, possibly dotted. There are no unenumerated continuous cross contacts. Interchanging the two chosen pieces is just a relabelling of an abstract quadrilateral, not an assumed symmetry of an ambient polytope.

The traversal terminates without a cutoff after 12,746 nodes. It gives **115 labelled union supports**, or **52 marked isomorphism classes**. Every deduplication has a saved full-pair bijection preserving the black nodes, opposite relation and labels. Finite tests use exact rational forest Schur elimination. Cycles are nonpositive already at their least allowed bonds; an edge of label 6 or 8 in a finite diagram must be an isolated rank-two component.

Disjoint orthogonal Lannér pieces are rejected because their spans supply two negative directions. Otherwise all extra union contacts are retained; no bridge-only assumption is made.

## Exact principal-subdiagram obstructions

The signature filter works in
\[
\mathbb Q\bigl(\sqrt3,\sqrt{2+\sqrt2}\bigr),
\]
using exact algebraic-field matrices and rational isolation of its real primitive element. For each finite-label core it eliminates a positive-definite corner and tests the remaining two-by-two Schur matrix. Both physical diagonals are negative. A positive determinant therefore proves two negative directions, independently of every omitted wall.

The other embedding sends \(\sqrt2\) to \(-\sqrt2\), and \(\sqrt{2+\sqrt2}\) to \(\sqrt{2-\sqrt2}\). Its Schur diagonals and determinant must be nonnegative. Fundamental-cycle parity also checks that cyclic products belong to \(\mathbb Q(\sqrt2)\); an odd occurrence of the independent label-6 or label-8 radicals violates this condition.

Of the 52 classes, **40 have two negative directions and one fails the cyclic-field condition**. The remaining classes are two dotted cores and nine finite-label cores. The latter have IDs
\[
5,9,14,17,32,37,41,47,50.
\]
These IDs refer to the explicit saved class list. Seven are singular; classes 14 and 37 are nonsingular. A principal core passing these tests is not yet a complete two-face star.

## Complete one-white contact domains

Any additional white node must have one of the finite contact patterns checked against the union core. For the nine finite-label cores, exact finite-corner and short-path tests leave **44** nonempty patterns. The field, physical-signature and conjugate-positive-semidefinite tests leave only **four**:

| Core | Permitted nonempty contact to one new white node |
|---|---|
| 5, 9, 17, 32, 41, 47, 50 | none |
| 14 | label 3 to black node 2 only |
| 37 | label 3 to black node 0 only; to black node 2 only; or to both |

All other contacts in each accepted row are right. For singular conjugate cores, every new row must annihilate the exact nullspace; this rejects 24 rows. The other rejections are saved with their exact field or Schur witnesses. The calculation imposes no bound on the number of further white nodes.

The seven singular cores are therefore **sealed**: every other white node is orthogonal to every node in such a core. Their affected finite corner components are consequently entire components, regardless of how many orthogonal white components exist elsewhere.

For the two dotted cores, the finite-corner calculation is independent of the dotted value. It leaves 12 and 6 patterns respectively. Apart from contacts to the two dotted endpoints, the only possibilities are a new 3-neighbor of black node 1, a new 3-neighbor of black node 3, or a common 3-neighbor of both. A single such neighbor produces a global affine fork with arm labels \(4,3,3\), forbidden by compactness. A common neighbor gives a \((3,3,4)\) Lannér triangle orthogonal to the dotted pair, hence two negative directions. Thus every external white node can contact the dotted core **only at its dotted endpoints**. This includes the common-neighbor case; it is not silently discarded as a pendant attachment.

## Arbitrarily large white supports: the finite-label cores

For core 14, a complete star cannot have no external contact: its nonsingular connected Gram would remain nonsingular after adjoining only orthogonal positive white blocks, whereas a complete quadrilateral star must have corank one. Its sole possible contact is a 3-neighbor at node 2. Adding that node gives a literal marked copy of sealed core 17. No further white node can contact this larger core. This closes every extension of class 14, including paths through whites initially having no contact with the old core.

For core 37, write \(a=0,y=2\). Its opposite contact is label 4, and both endpoints touch the middle of the fixed C5 by label 3. Every external contact is at \(a\) and/or \(y\), by the table. The white graph is a forest. An external connected component cannot have two contacts to the same endpoint: that would create a cycle in a finite corner. Two distinct components at the same endpoint would make a second branching vertex alongside the existing D5 branch, also impossible in a finite diagram.

In a corner containing \(a\), an attached component extends the short arm of D5. A finite connected diagram containing this branch is simply laced of type D or E, has only one branching vertex, and permits only a simple path of length \(0,1,2,3\) at that arm. These yield D5, E6, E7, E8. Length four contains affine E8. The same statement holds at \(y\). This classification is the all-size attachment argument; the four finite cases and first affine obstruction are exact checks, not a size cutoff.

If an external component connects \(a\) to \(y\), its white part is therefore a path with one, two or three nodes. Together with the opposite 4-bond it gives a C3, C4 or C5. The existing C5 and all extra old contacts remain present. The resulting full union is respectively a saved literal copy of sealed core **41, 47 or 50**. The checker records the complete marked pair bijections. Hence these extensions also stop without further attachments.

If no component connects \(a\) and \(y\), only separate pendant paths remain. Eliminating them changes the two unit diagonals to
\[
p,q\in\left\{1,\frac34,\frac23,\frac58\right\}.
\]
For the original norm-two opposite contact \(r=\sqrt2\), direct full-star rank gives
\[
(p+q-4pq)+(1-p-q)\sqrt2=0.
\]
Both coefficients must vanish, forcing \(p+q=1\), while the allowed tails give \(p+q\ge5/4\). Thus a nonsingular class-37 core cannot be completed only by such pendants. This covers all its external components, including branches and paths through vertices not initially adjacent to the core.

## Arbitrarily large white supports: the dotted cores

Both remaining cores have a C5 \((b,u,v,w,d,b)\) with its only 4-bond on \(b,d\), and a dotted pair \(a,y\). Node \(a\) touches \(v\) by a 3-bond; node \(y\) is either right to \(v\) or also touches it by 3. The complete one-white calculation proved that no external contact is possible except at \(a,y\).

If no external white component connects both endpoints, every component is pendant at one endpoint or orthogonal. The [independently reviewed pendant classification](fresh_sqrt2_dotted_tail_20261015.md) then applies without an attachment cutoff. It leaves only the untailed single-contact dotted star, its one-tail label-8 version, and the untailed double-contact dotted star. Their affected complete components are D5 or E6. The label-8 version is already sealed core 17.

It remains to exclude a white component connecting both dotted endpoints. Eliminate **all** external white nodes. Their matrix is positive definite with nonpositive off-diagonal entries, so its inverse is entrywise nonnegative. In unit normalization the resulting endpoint diagonals satisfy \(0<p,q\le3/4\): each endpoint had at least one nonright contact of cosine at least \(1/2\), and the corresponding diagonal inverse is at least one. The off-diagonal correction is strictly negative because the two endpoint contact supports meet the same connected white component. Thus the effective norm-two opposite magnitude satisfies
\[
r_{\rm eff}>r>2.
\]
The endpoint bonds of the connecting component need not all be 3; labels 4, 6 or 8 at an endpoint only strengthen these inequalities.

The two exact rank identities from the pendant package remain valid for these reduced diagonals and increased effective contact, because elimination creates no contacts elsewhere in the core. In the single-contact case,
\[
r_{\rm eff}^2=q(4p-1+\sqrt2)
\le\frac34(2+\sqrt2)<4,
\]
contradicting \(r_{\rm eff}>2\).

For the double-contact case put \(t=(1-\sqrt2)/2<0\). The positive root is
\[
r_{\rm eff}=-t+\sqrt{(2p-t)(2q-t)}.
\]
It increases with both positive variables. Therefore
\[
r_{\rm eff}\le\frac12+\sqrt2<2,
\]
again a contradiction. The other root is negative. This excludes all multiply attached exterior components and justifies reduction to the reviewed pendant theorem.

We have now exhausted every compact quadrilateral with no distance-two marked corner. No assumption that the white support was already contained in the initial minimal-pair union remains.

## Face weights and the dimension consequence

Assign weight 1 to every marked pair at graph distance at most two in its finite vertex diagram. Add the following nonnegative bonuses, **only when the indicated connected component is entire**:

| Entire finite component | Additional marked pairs and weights | Total bonus in that component |
|---|---|---:|
| F4 | the two ends of the path, weight 1 | 1 |
| D5, arms 1,1,2 | long-arm tip with either short-arm tip, weight 1/2 each | 1 |
| E6, arms 1,2,2 | near node on one length-two arm with far node on the other, weight 1/2 each | 1 |
| E7, arms 1,2,3 | near node on the length-three arm with the length-two-arm tip, weight 1/4 | 1/4 |
| E8, arms 1,2,4 | near node on the length-four arm with the length-two-arm tip, weight 1/4 | 1/4 |

The exact corner check gives the following payments for the sealed cores: class 5 has one F4 corner; class 9 has two; class 32 has four; class 17 has two E6 corners; class 41 has four; classes 47 and 50 have four E7/E8 corners respectively. The dotted stars have two or four D5 corners. Thus every exceptional quadrilateral receives at least one unit. Other quadrilaterals already received one from a short marked pair.

Every triangular face receives at least two units from the basic distance-two rule. A minimal Lannér subset of its full star contains all three black nodes. If it has order three, its connected graph has at least two nonright edges. For C4/C5 the three cyclic gaps between the black nodes are positive integers summing to four or five; at least two gaps have length at most two. Their paths have only white internal nodes, so each occurs in the corresponding finite corner. Larger faces receive nonnegative weight, which is at least \(5-k\) for \(k\ge5\). Therefore every two-face with \(k\) sides has total weight at least \(5-k\).

For an irreducible finite component of rank \(r\), the basic count is \(2r-3\) on a path and \(2r-2\) on D/E trees. After bonuses, D5 and E6 have deficit \(2r-w=1\); E7/E8 have deficit \(7/4\); all other components have deficit at least two (A1 also has deficit two). For total rank \(n\ge9\), a single component has deficit at least two, because all exceptional components have rank at most eight; two or more components have total deficit at least two. Consequently every vertex receives at most
\[
2n-2.
\]

The standard h-vector/face-incidence inequality used in the [reviewed compact-bound package](fresh_compact_dimension_bound_20261015.md) makes the average vertex weight strictly greater than
\[
\begin{cases}
n(n-6)/8,&n\text{ even},\\
(n-1)(n-5)/8,&n\text{ odd}.
\end{cases}
\]
For odd \(n\ge21\), the lower expression minus \(2n-2\) is \((n-1)(n-21)/8\ge0\). For even \(n\ge22\), it is \((n(n-22)+16)/8>0\). Strictness also excludes equality at \(n=21\). This gives the stated dimension-at-most-20 consequence. The separate [weight certificate](fresh_sqrt2_weight_bound_20261015.md) and its [independent review](fresh_sqrt2_weight_bound_20261015_review.md) also check the complete component charges and scalar inequality.

## Verification scope

All four primary generation commands and their `--verify` replays pass. Independent generation and saved replay also pass. The finite signature and extension tests use exact number-field arithmetic, not numerical eigensolver tolerances. The source manifests bind the complete Lannér catalogue and the independently reviewed dotted-tail theorem. The unbounded white-component argument and the h-vector implication are mathematical proofs, not extrapolations from the finite graph tables. The immutable closure JSON retains its historical pre-review status string; the linked, separately hash-bound review records the completed independent validation.

No earlier frozen source or certificate is modified. Independent hostile review of the complete support enumeration, attachment closure and final scope reports no gap. This package is now frozen.
