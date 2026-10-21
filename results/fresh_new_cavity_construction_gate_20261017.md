# Exact mandatory-endpoint obstruction for a four-center pentagon disk

Exact generation and a fresh saved-output replay both passed the complete finite-label domain of the specified four-cap disk. This is an obstruction to that topology with its stated retained walls, not to arbitrary compact or finite-volume completion. No new chamber is promoted.

Start from the frozen [63-wall, 199-finite-vertex control](fresh_nonarithmetic_release49_62_replacement_20261017.md) and its [complete S6 boundary](fresh_nonarithmetic_release49_62_boundary_20261017.md). Here S6 is the positive-definite H4+A2 support `[0,1,2,3,5,6]`. Delete exactly `[12,47,52,64,65]`. The changed exterior pentagon is `[7,8,44,63,30]`; the retained chord 30–44 and finite triangle 30–44–63 lie outside the proposed disk. The inherited deletion ledger removes eleven old finite-vertex certificates and six certified compact two-faces. Its exact supports are saved, with no assumption that a replacement restores them.

The new cap arcs are A over 7–8, R over 8–44–63, B over 63–30 and C over 30–7. The interior cycle A–R–B–C–A is triangulated by either A–B or R–C. Both versions necessarily contain the finite corner S6+8+A+R. The proof excludes the A,R pair, so it does not depend on the remaining labels or normal equations for B,C. The [numerical predecessor](fresh_new_cavity_construction_20261017.md) motivated this necessary-corner check but is not used to certify a sign or branch.

Every actual finite endpoint at the old edge 7–8 contains the H4 component on A2+7+8. No larger finite irreducible Coxeter diagram contains H4. Thus A is right to all of S6, 7 and 8. For R, the two required finite corners S6+8+44+R and S6+44+63+R first force rightness to H4 and one of the seven finite A2 attachment profiles. The complete finite classification of the remaining rank-five diagrams gives 96 recipes whose three old labels lie in 2 through 5, and ten high-label families. A label at least 6 belongs to an isolated I2 component in each required finite corner; the symbol 6 therefore represents every such label, including arbitrarily large ones. Applying the same classification at S6+8+A+R gives 285 cases: 271 low A–R labels and fourteen high-label families. The small graph tests also explicitly verify that A's only profile at 7–8 is right to A2, 7 and 8.

The exact geometry is computed in the Lorentzian four-dimensional complement of S6. If a normalized normal has A2 contact vector c, its projected square is

\[
q(c)=2-c^{\mathsf T}\begin{pmatrix}2&-1\\-1&2\end{pmatrix}^{-1}c.
\]

All other native products are recovered by adding the corresponding A2 projection product. The source reconstructs the literal projected Gram of `[8,44,63]` and its inverse. Projecting old wall 30 off this frame gives a complementary vector K with strictly positive square. For every contact vector t the complete R norm equation is

\[
R_{\pm}=\sum_i (G^{-1}t)_i n_i
 \;\pm\sqrt{\frac{q(c)-t^{\mathsf T}G^{-1}t}{K^2}}\,K.
\]

Both signs are retained. A negative remaining square rules out a real branch because the complementary line is positive. Intervals crossing zero are conservatively intersected with the nonnegative half-line before taking the square root; this does not discard a real root.

For each R branch not already contradicted by a retained native product, the second frame is `[7,8,R]`, with A contacts `(0,0,-2cos(pi/m))` and square 2. Its inverse and a strictly positive complementary square are certified before interpreting its norm equation. In a Lorentz space of index one this proves that the nondegenerate three-frame is Lorentzian. The same two-branch formula therefore covers every possible A. A frame for which these interval checks are inconclusive is subdivided; it is never accepted as an obstruction merely because inversion or a sign test failed.

Each high cosine is initially relaxed to the closed interval `[sqrt(3),2]`, including the parallel limit as an enlargement. Exact rational bisections cover the entire initial box. A leaf is certified by either a strict obtuse retained-wall product or a strict gap between consecutive Coxeter products. Such a gap may occur in a prescribed high contact itself, or in a native product of R or A with an old retained wall. Consequently the theorem covers all Coxeter labels; it does not claim an obstruction for every real contact parameter. The source searches consecutive thresholds through label 65 to locate witnesses. This search is not a candidate-label cutoff: a branch without a strict certificate remains unresolved and makes generation fail. Every displayed witness is a strict exact inequality, so equality and boundary values are never silently rejected.

The complete output has 2,305 tested interval boxes and 1,295 certified leaves, with maximum subdivision depth 21. There are 452 leaves excluded by a prescribed-angle gap and 843 leaves checked through the norm frames. The latter certify 1,089 R-branch native obstructions and 597 surviving R branches for which both A branches are excluded. Thus all 1,194 corresponding A branches are rejected. No unresolved leaf remains.

The [source](../scripts/fresh_new_cavity_construction_gate_20261017.py) and [certificate](fresh_new_cavity_construction_gate_20261017.json) bind the frozen physical geometry, finite-diagram classification and predecessor inputs. They reuse the established outward rational interval engine; this package is not an independent implementation of that engine. Root independently read the nested-frame argument, finite-profile reduction, positive-complement guard and gap scope. That review is mathematical/source review, not a second interval-arithmetic replay or machine formalization.

Generation handle 34082 and fresh replay handle 1977 both exited zero. Syntax, all eleven direct source hashes and all five local report links passed. The numerical predecessor independently passed its own fresh saved-output replay.

This closes one new four-center incidence template and both of its central diagonals. Other four-cap assignments, a center with no old boundary edge, additional intersections with old sides, further releases or deformations, and walls missing S6 remain outside the claim. The exact 63-wall control, its strict escape and its infinite-volume status are unchanged.
