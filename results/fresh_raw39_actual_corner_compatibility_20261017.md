# Actual-corner compatibility within the raw39 square family

The new family gives four violated necessary columns at three vertex forests. Repairing only those vertex allowances lowers the [new local-and-joint witness](fresh_mixed_actual_corner_full_model_20261017_joint_ports.md) from12.38638766981714 to **12.382070564084598**, still above12. All4,661 face masses remain fixed. This is an exact fixed-witness repair, not an optimized value, a target12 feasible rule or a dimension theorem.

The [source](../scripts/fresh_raw39_actual_corner_compatibility_20261017.py) and [certificate](fresh_raw39_actual_corner_compatibility_20261017.json) cover all31 positive classes of canonical core29, whose literal raw source is39. One class uses the existing explicit provenance map; the other30 have exactly one saved source context each. The complete source-context list is checked before lifting a merged row. No arbitrary realization is chosen, and there are no omitted ambiguous positive raw39 classes in this input.

The core is a seven-root local square with four finite corners and dotted normalized square

\[
\tau=(41+17\sqrt5)/2.
\]

An orthogonal rank8 finite spectator forest raises each corner to rank13. Each corner retains five raw roots, with finite component types A1+A4, A5, A5 or2A1+A3. The source reconstructs all124 actual corner occurrences,88 resulting whole vertex forests and1,246 distinct literal orientations. Every complete active component embeds into a whole matching vertex component, over every identical component copy and every internal weighted automorphism. Passive components have zero external covectors and may be quotiented. A full literal star is identified by its unordered pair of complete external contact covectors in the common labelled finite-vertex frame. Repeated keys must agree in occupied pair and native omitted product. Complete orbit keys are disjoint; a corner orbit belongs to one actual Q-face class. Coefficients count actual corner occurrences, not the number of orientations.

At a finite vertex, let G be its positive definite incident norm2 Gram matrix. For an external actual wall with contact covector p, its timelike height is the positive number

\[
h_p=p^TG^{-1}p-2>0.
\]

The future outward sign is fixed by the common vertex. Two external actual wall normals therefore have forced product

\[
g(p,q)=p^TG^{-1}q-\sqrt{h_ph_q}.
\]

The literal native product for each reconstructed square is verified by this formula. Identical covectors represent the same wall and are skipped in incompatibility tests. A positive cross product is impossible for distinct walls of an acute Coxeter chamber, as is a product strictly between−1 and0 for norm2 walls. Distinct complete square stars also cannot occupy the same actual two-port pair at one simple vertex: that pair fixes a unique two-face and its full wall star.

The arithmetic tests explicitly use the fixed ambient ground field K=Q(sqrt5), not arithmeticity over an unspecified field. Here tau is nonsquare in K, and every cross product has form

\[
g=a+b\sqrt\tau-\sqrt h,\qquad a,b,h\in K,\quad h>0.
\]

The source first determines the exact physical sign, then decides whether g² belongs to K. If sqrt(h)=r lies in K, write A=a−r,B=b. If instead sqrt(h/tau)=r lies in K, write A=a,B=b−r. In either quadratic case, g² is in K exactly when AB=0. If both square-root tests fail, K(sqrt(tau),sqrt(h)) is biquadratic; expansion in its independent basis shows g² lies in K exactly when a=b=0. These exhaust the cancellation cases. A surviving square must be an algebraic integer in O_K by the native two-cycle criterion. For a finite acute cross product, the remaining fixed-K angle alphabet is also checked. These are necessary exclusions; no converse realization assertion is used.

All21,115 unordered literal-star pairs are tested. Exactly21,047 are excluded and68 remain conservatively compatible. The recorded first reasons are2,128 shared square-port pairs,3,197 positive distinct-wall products,9,626 nonintegral two-cycles and6,096 squares outside K. Every one of the124 corner orbits is internally incompatible, hence has capacity1. There are nine universal incompatibility edges between such orbits. Maximal orbit cliques and verified proper colorings of the conservative compatibility graphs give129 exported columns. Unknown pairs remain edges in that compatibility graph. Simultaneous actual stars form a clique, so a proper coloring is a valid upper capacity even if it is not optimal.

The four violations reduce to three allowance increases:

| Whole vertex forest | New required allowance | Previous allowance |
|---|---:|---:|
| 6A1+A3+D4 | 579797/10000000000 | 2931717/100000000000 |
| 5A1+A3+D5 | 647797/10000000000 | 143887/2500000000 |
| 2A5+H3 | 12822341/5000000000 | 1124143/500000000 |

The last forest has both a capacity-one single-orbit column and a redundant capacity-two union column. After taking per-vertex maxima, all129 new loads pass exactly. Only vertex allowances increase, so the preceding complete local simplex covers, all17 selected joint families, every inherited cap and incidence inequality, and all432,021 already checked model columns remain valid by monotonicity.

The exact final ratio is

\[
\frac{247898037929594052000000000000000}
{20020725665111814911449997327399}
=12.382070564084598\ldots>12.
\]

This package is confined to within-family raw39 compatibility at the stated positive support. Cross-family raw39/raw23/raw25 interactions, other Q classes, missing face and forest domains and full ambient compatibility remain outside it. No optimizer is invoked.

Generation89935 and fresh saved replay9555 both terminated exit0 with `PASS_EXACT_RAW39_UNIQUE_PROVENANCE_FINITE_VERTEX_CORNER_CAPACITIES`. Replay reconstructs the actual corners, orientation domains, pair certificates, capacities and exact repair without an optimizer. The three-branch field-membership argument is the same previously reviewed algebra used by the frozen [mixed raw23/raw25 package](fresh_mixed_raw23_raw25_corners_20261017.md); that prior review is not a new raw39 projection or census replay.
