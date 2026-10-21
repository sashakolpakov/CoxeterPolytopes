# Mixed raw23/raw25 actual-corner compatibility

The exact fixed-witness repair lowers the saved bound from **12.38203170968089 to 11.863463132596562**. All 4,932 face masses stay fixed. The 65 increased vertex allowances are forced by genuinely mixed raw23/raw25 columns; the columns involving only one family require no increase. This removes this particular obstruction above 12. It is neither an optimized bound nor a feasible weighting rule, and it gives no new dimension theorem.

Primary [source](../scripts/fresh_mixed_raw23_raw25_corners_20261017.py) and [saved certificate](fresh_mixed_raw23_raw25_corners_20261017.json). Generation and a fresh `--verify` reconstruction both completed with exit 0 and `PASS_EXACT_MIXED_RAW23_RAW25_FINITE_VERTEX_CORNER_CAPACITIES`. No optimizer was invoked. The input and mathematical scope remain compact arithmetic with fixed ground field K = Q(√5).

## Complete represented domain

The baseline is the frozen [raw25 certificate](fresh_raw25_actual_corner_compatibility_20261017.md). Its 4,932 face masses and all vertex allowances are read exactly. The new source reconstructs both cohorts from their literal core and spectator records:

| Family | Positive face classes | Actual prism corners | Literal orientations | Whole vertex forests |
|---|---:|---:|---:|---:|
| raw23, canonical core 19 | 57 | 570 | 6,422 | 308 |
| raw25, canonical core 21 | 55 | 440 | 6,260 | 296 |
| Union | 112 | 1,010 | 12,682 | 379 |

There are 225 common forests and exactly 227,064 cross-family orientation pairs. The current raw23 cohort removes zero-mass rows 26457 and 26465 and adds row 26259 (provenance 440, context 27280, eight A1 spectators) and row 26372 (provenance 472, context 27393, two A1 and three I2(5) spectators).

All 57 raw23 classes and 19 raw25 classes use the inherited injective actual-provenance map. The other 36 raw25 rows have exactly one saved source context each, whose canonical core and spectator rank are verified. No representative is chosen arbitrarily for a row with several saved provenances. This is coverage of the stated represented classes; it does not assert coverage of every possible arithmetic face.

For each corner, every complete active finite component embeds into every matching whole component of the canonical rank-13 vertex forest, with every internal weighted automorphism and choice among identical spectator copies. Passive components have zero external contacts. The unordered pair of complete external contact covectors is the literal full-star key. Identical keys are merged, and their occupied square pairs and prescribed omitted-root product must agree. Different complete orbit keys are explicitly checked to be disjoint. Each orbit has one face class; the global coefficient is the number of actual corners in that orbit, not the number of its orientations. Raw23 has ten actual prism vertices and raw25 eight; the quad partitions determine these vertices, rather than an unrestricted list of spherical deletions.

## Exact common-frame geometry and field gate

Let v be the future unit timelike vector of a finite simple vertex, let G be its positive-definite normalized incident Gram matrix, and let p be the complete covector of an external norm-2 wall. The wall has the unique outward expression

\[
r_p=s_p+\sqrt{p^T G^{-1}p-2}\,v,
\qquad s_p=\sum_i(G^{-1}p)_i n_i.
\]

The square root is strictly positive because the wall is external. Thus two external walls in the same complete vertex frame have forced product

\[
(r_p,r_q)=p^T G^{-1}q-
\sqrt{(p^T G^{-1}p-2)(q^T G^{-1}q-2)}.
\]

The source checks positive heights and recovers the prescribed native product within every literal star. Identical external covectors denote the same actual wall and are skipped in incompatibility tests.

Put τ = (7+√5)/2. The raw23 dotted contact lies in K√τ; the raw25 dotted square is 6+2√5 = (1+√5)², so its contact lies in K. The common projection therefore has the form a+b√τ with a,b in K, while the height product h is positive in K. All signs of g = a+b√τ−√h are decided exactly, with the sign checked before squaring.

The native two-cycle g² must belong to K and be integral in O_K. Membership is decided completely as follows, with positive physical square roots:

* If √h = r belongs to K, write g = (a−r)+b√τ.
* Otherwise, if √(h/τ) = r belongs to K, write g = a+(b−r)√τ.
* In either case g = A+B√τ has square in K exactly when AB = 0, since τ is not a square in K.
* If neither square root belongs to K, K(√τ,√h)/K is biquadratic. The non-K coefficients in g² are 2ab, −2a, −2b on √τ, √h and √(τh). Hence g² belongs to K exactly when a=b=0, in which case g²=h.

Integrality is tested only after membership is established. For a finite acute contact, the fixed-K squared-product alphabet is {0,1,2,(3+√5)/2,3,(5+√5)/2}, corresponding to labels 2,3,4,5,6,10. The fixed-field premise and its distinction from integrality alone are inherited from the [ground-field gate](fresh_sqrt5_spectator_groundfield_gate_20261017.md). This calculation does not apply the raw25-only K-projection shortcut to mixed products.

Two distinct complete stars are also incompatible if they occupy the same actual square port pair: that port pair determines one actual two-face and its complete full wall star. Exact positive products and products strictly between −1 and 0 exclude pairs independently of the arithmetic tests.

## Pair census and valid capacities

All 535,065 unordered literal pairs in the 379 forests are checked. The certificate excludes 516,727 and conservatively retains 18,338 as compatible. Of the cross-family pairs, 216,446 are excluded and 10,618 remain allowed. Exclusions use:

| First certified reason | Pairs |
|---|---:|
| Shared actual square port pair | 85,419 |
| Positive external-wall product | 190,648 |
| External-wall product in (−1,0) | 132,808 |
| Nonintegral native two-cycle | 101,544 |
| Native two-cycle outside fixed K | 6,308 |

There are 600 internally incompatible corner orbits and 435 universally incompatible pairs of such orbits. Their cliques supply capacity-one columns. Further columns use explicit proper colorings of the conservative compatibility graph on every single orbit and on the union of all represented orbits in each forest. Actual simultaneous stars form a clique in that graph, so the number of colors is a rigorous capacity. Every color class is checked against the exact pair certificates. No optimal coloring or complete geometric classification is asserted; unexcluded pairs remain compatible.

The resulting 1,251 columns include 245 genuinely mixed columns. Exactly 67 columns violate the baseline, all mixed. They require 65 vertex-mass increases. For example, at forest A1+A4+A5+H3 (vertex 3165), one exported capacity is

\[
\alpha_{26410}+2\alpha_{26477}+\alpha_{26947}\le\mu_{3165}.
\]

At A1+A2+A5+H3 (vertex 2622), another is

\[
\alpha_{26392}+\alpha_{26413}+\alpha_{26477}+\alpha_{26947}\le\mu_{2622}.
\]

The multiplicity two counts two actual corners of the relevant face, not two presentations at one vertex. The saved literal orbits and pair witnesses certify the corresponding capacity-one cliques.

## Exact repair and preserved scope

For every exported column with positive cost c at vertex forest v, set μ_v to the maximum of its previous value and the exact load divided by c. All α remain unchanged. Every new column is then recounted exactly. Only nonnegative vertex allowances increase, so all inherited model, signed-transfer, signed vertex-cost incidence, complete fixed-mass local-packing, and previously checked joint-family inequalities remain valid by monotonicity. No changed-mass cover is imported.

The unchanged objective is

\[
N=4132201186871/330000000000.
\]

The repaired normalizer and ratio are

\[
D=5705370740772271626674000656033/
5405400000000000000000000000000,
\]

\[
N/D=
\frac{67685455440946980000000000000000}
{5705370740772271626674000656033}
=11.863463132596562\ldots.
\]

The single-family-only normalizer is exactly the baseline normalizer, so the threshold crossing comes entirely from mixed compatibility. The certificate contains all original face masses, baseline and repaired vertex masses, exact per-column loads, explicit colorings, orbit cliques, native-product checks and source hashes. Its compact pair schema stores one certified reason for each excluded pair and shares repeated radical expressions within each forest.

The independent narrow algebraic review confirms the complete square-class membership split. The root's separate full-source and mathematical read also passes the common-frame maps, unique provenance, same-wall treatment, actual-corner multiplicities, clique/coloring capacities and monotone repair. These are scoped mathematical/source reviews, not independent orbit or projection replays, and are distinct from the producer's full generation and fresh reconstruction. A new finite optimization would need to retain the earlier exported matrix and these positive-cost columns. Even an exact target-12 primal for that finite model would still require the missing global face, field and weighting coverage before any dimension conclusion.
