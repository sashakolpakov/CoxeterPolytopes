# Independent review of the exceptional distance9 supporting extensions

The proposed residual-space argument closes the supporting-white extension problem on the supplied69 exceptional contexts. Independent exact reconstruction gives nine full-core types,24 relevant-root cases and46 contexts with no relevant core-contact node. The resulting marked-component order is at most16. No gap was found in the indirect-tail or no-bridge argument under the stated hypotheses.

The [independent source](../scripts/fresh_distance9_supporting_extensions_review_20261017.py) imports no producer arithmetic. Its [certificate](fresh_distance9_supporting_extensions_review_20261017.json) reconstructs normalized Grams directly from the frozen [quasi-Lanner census](fresh_high_rank_quasilanner_census_20261017.json), using SymPy rational matrix arithmetic. The69 literal contexts and selected corners are explicit inputs in the source. This review checks those contexts, rather than independently proving that they exhaust all high-core overlaps. The whole distance9 classification, unknown core cross-edge exclusions, low-order cases and dimension composition remain separate obligations.

Each reconstructed core H has inertia(order(H)−1,1,0). Consequently its orthogonal complement in the ambient Lorentz space is positive definite. For every selected corner, the source recomputes its connected finite or affine blocks, their exact inverse diagonals or positive kernels, the complete necessary core-contact set, and the marked component. Affine blocks permit no external contact: their positive null relation and acuteness force all such products to vanish. Finite blocks give the usual inverse-diagonal bound2. The selected corner places every feasible projection-one contact outside the marked component.

Write J for the inverse normalized full-core Gram. Suppose the short-path alternative has not already paid the face. The saved end partition has maximum left distance plus maximum right distance plus2 at most9. A supporting white root contacting both ends therefore supplies a short path. Within either end, all entries of J are positive, and the exact minimum projection for any two distinct contacts is greater than2. Since every nonzero acute Coxeter contact has magnitude at least1, no multicontact root remains. Thus every core-contact root has a single contact. The candidate nodes satisfy1≤J_ii≤2. If J_ii>1, only a simple bond is possible; if J_ii=1, a simple or doubled bond is possible. Larger labels have projection square greater than2.

For a simple-contact root at node i, write its core projection as v_i and its positive-space residual as g_i. Then

\[
v_i^2=J_{ii}=p_i,\qquad g_i^2=2-p_i=\rho_i,
\qquad r_i\cdot r_j=J_{ij}+g_i\cdot g_j.
\]

Cauchy bounds the residual product. The exact pair tests either force a positive actual product, force a nonzero native edge giving a path of length at most9, or place the actual product strictly above−1. In the last case the discrete Coxeter gap and acuteness force the product to be0. All doubled-contact candidates at projection-one nodes are separately checked: their residual vanishes, so the sign of the exact core cross product either excludes coexistence or forces a short path. These checks also exclude simultaneous relevant roots, apart from the single relevant root being analysed.

Let G be that relevant root, with0≤rho<1. There is at most one other allowable core-contact type C when G is present, and C contacts a node outside the marked component. The actual roots G,C must be orthogonal. If its simple-contact projection is p_C>1, two distinct roots at the same contact node would have actual product at least2p_C−2>0. If p_C=1, two such roots could be acute only by being orthogonal with opposite unit residuals. Both cannot have the same strictly negative residual product with G. Hence at most one C occurs. Identical complete roots are a single actual wall, not two occurrences.

The four exact residual possibilities are:

| p_G | p_C | projection cross product P | residual determinant | inverse GG entry |
|---:|---:|---:|---:|---:|
| 4/3 | 4/3 | 2/3 | 0 | — |
| 3/2 | 1 | 1/2 | 1/4 | 4 |
| 7/4 | 1 | 1/2 | 0 | — |
| 5/4 | 1 | 1/2 | 1/2 | 2 |

In each case the residual Gram is S=[[rho,−P],[−P,rho_C]]. If det(S)=0, its kernel has strictly positive coordinates, so a positive combination of the two residual vectors is zero. Any white root orthogonal to the full core has nonpositive products with both G and C; the positive relation forces both products to be zero. No indirect attachment can meet either root.

If det(S)>0, S is a positive definite Stieltjes matrix and S^{-1} has nonnegative entries. A norm2 white root E orthogonal to the core, but joined to G, has contact magnitude at least1. Its projection on span(g,c) therefore has squared norm at least(S^{-1})_GG. If that entry is greater than2, E is impossible. If it equals2, every inequality must be equality: E has contacts(−1,0), no component outside span(g,c), and

\[
E=-2g-(S^{-1})_{CG}c.
\]

For any other core-orthogonal white Z, acuteness gives g·Z≤0 and c·Z≤0, hence E·Z≥0. Acuteness again forces E·Z=0. Also E·C=0. The equality determines E uniquely. Thus the one permitted branch is sealed; it cannot bridge to another outside component or continue through further core-orthogonal whites.

If C is absent, all other supporting whites coexisting with G are orthogonal to H. Let E be a whole connected component among them, of rank m, and let b be its nonnegative contact magnitudes with G. Positive-space projection gives b^T G_E^{-1}b≤rho. For every connected positive definite norm2 Coxeter Gram of order m,

\[
(G_E^{-1})_{ii}\ge m/(m+1).
\]

This bound does not assume a bounded rank classification. Choose a simple spanning tree. Its nonnegative adjacency matrix is entrywise at most the original one, so positive definiteness and the Neumann-series inverse comparison reduce to that tree. Root the tree at i. The Schur complement gives inverse diagonal1/(2−sum of rooted child inverse diagonals). By induction, child subtrees of sizes m_j contribute at least m_j/(m_j+1), whose sum is at least(m−1)/m. This yields m/(m+1). The base case is1/2. The same inverse comparison handles larger bond magnitudes.

A contacted component consequently satisfies m≤floor(rho/(1−rho)). Two distinct contacted components would each cost at least1/2, contradicting rho<1. The bound controls the entire indirect component, including roots having no direct contact with G. It is not merely a bound on the first attached node. In the present cases the maximum alone-tail rank is3.

The no-bridge conclusion now follows in both alternatives. With C present, no core-orthogonal route can continue from G except the sealed single root just described; G and C themselves are orthogonal. With C absent, at most one complete finite component of the stated rank can attach to G, and all other core-contact root types have already been excluded or would force a short path. If there is no relevant core-contact node, nothing can enter the marked core component at all. The source combines these bounds with each actual preferred component and obtains order≤16 in all69 contexts.

Generation and a fresh saved replay both pass `PASS_INDEPENDENT_DISTANCE9_SUPPORTING_EXTENSION_CONTEXTS`. The replay repeats the exact69-context construction and all residual tests. This is an independent local arithmetic and mathematical review; it makes no independent complete-overlap enumeration or new dimension claim.
