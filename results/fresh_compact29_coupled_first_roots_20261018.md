# A fixed-witness incompatibility at the E7³+E8 vertex

The frozen individual small-face witnesses cannot all coexist at one \(E_7^3+E_8\) finite vertex. Five specified E7 recipes, when placed on two distinct E7 factors, force a mutual product strictly between Coxeter labels11 and12. This gives75 excluded literal pair placements in a complete domain of3,087 placements of the **first added roots** of the49 saved E7/E8 recipes. It is a new compatibility test on those fixed witnesses, not an exclusion of the vertex type, a complete classification of its face stars, or a compact-dimension improvement.

The [producer](../scripts/fresh_compact29_coupled_first_roots_20261018.py), [certificate](fresh_compact29_coupled_first_roots_20261018.json), and [standalone checker](../scripts/check_compact29_coupled_first_roots_saved_20261018.py) bind the [exact79 individual witnesses](fresh_compact29_small_faces_20261015.md) and the [rank29 savings target](fresh_compact29_savings_frontier_20261018.md). The original weights and recipes remain frozen. No arithmeticity hypothesis is imposed.

## Forced products in one common Lorentz space

Let \(T=E_7\oplus E_7\oplus E_7\oplus E_8\) be the positive-definite span of the29 old unit-normalized roots (each root has square2). Its orthogonal complement in the ambient signature\((29,1)\) space is a timelike line. Choose its unit future vector\(e\), with \(e^2=-1\); the old vertex lies on the ray\(e\).

For a first added root\(c\), its products\(b\) with the old roots determine its projection\(p\) onto\(T\), with \(p^2=b^\mathsf{T}G_T^{-1}b\). Norm2 forces

\[
c=p+h e,\qquad h^2=b^\mathsf{T}G_T^{-1}b-2.
\]

Each saved recipe has\(h^2>0\). Its old vertex must remain strictly inside the new wall, so\(c\cdot e=-h<0\): the positive sign of\(h\) is mandatory. The opposite sheet would exclude the very vertex whose face star is being extended.

For first roots\(c_i,c_j\) supported on distinct old irreducible factors, their old projections are orthogonal. Their mutual product is therefore forced, without a free label:

\[
c_i\cdot c_j=-h_i h_j,\qquad (c_i\cdot c_j)^2=h_i^2h_j^2.
\]

This argument uses the actual common old vertex and ambient dimension. It is not obtained by assigning unrelated negative directions to separately realizable stars. Extending each saved recipe orthogonally over the other three original factors is its explicit declared embedding.

## Complete declared domain and the native gap

There are21 E7 recipes and28 E8 recipes, one for each marked pair in the frozen source. The first added root is the triangle's\(c\), or the quadrilateral's\(c\); the latter includes its saved\(a\)-\(c\) contact. All49 old-product vectors are reconstructed exactly in\(\mathbb Q(\sqrt2,\sqrt5)\), with all positive height-squares saved.

The three E7 component slots contribute\(3\cdot21^2=1323\) literal pairs. Their symmetry quotient has\(\binom{22}{2}=231\) rows. The three E7/E8 slot pairs contribute\(3\cdot21\cdot28=1764\) literal pairs, with588 quotient rows. Thus all819 quotient rows represent exactly3,087 literal pairs. There is only one E8 factor, so an E8/E8 pair is not part of this cross-factor domain. Same-factor pairs and second quadrilateral roots are not tested here.

The five E7 marked pairs are

\[
(2,4),(2,5),(3,5),(4,5),(5,6).
\]

Their saved first roots all have

\[
h^2=\frac{1+3\sqrt5}{4},\qquad
h^4=\frac{23+3\sqrt5}{8}.
\]

The first four recipes are triangles; the last is a quadrilateral. Two such recipes on distinct E7 factors would require the squared mutual product\((23+3\sqrt5)/8\). Directed rational Machin-series and cosine-series bounds prove

\[
4\cos^2(\pi/11)<\frac{37}{10}
<\frac{23+3\sqrt5}{8}<\frac{93}{25}
<4\cos^2(\pi/12).
\]

For distinct norm2 Coxeter roots the squared magnitude of a nonright finite product is\(4\cos^2(\pi/m)\), strictly increasing with the integer\(m\ge3\). The displayed strict gap excludes every finite label, including arbitrarily large labels; the value is also below4, so it is neither parallel nor dotted. This is an exact native-product obstruction. It uses no field classification or numerical label selection.

These five-by-five choices give25 exclusions for each of three pairs of E7 slots:75 literal exclusions, represented by15 quotient rows. Every other first-root pair has squared product strictly greater than4, hence a negative strictly dotted product. The3,012 dotted literal pairs merely pass this one necessary gate. No claim that their complete face stars coexist follows.

## Consequence and remaining search

The independent local witnesses from the earlier packet cannot be assembled wholesale at one finite vertex. In particular, a compact-polytope argument cannot choose those first-root recipes independently on all three E7 factors. This is a concrete coupled-incidence obstruction that the individual finite-corner/inertia tests do not detect.

Alternate recipes may change the heights or contact supports. The present source does not classify them, and it does not prove that either selected angle has to use its frozen recipe. Even among retained first-root pairs, the remaining roots, mutual products, full small-face persistence, common interior, and global completion remain to be checked. Consequently there is no new face-charge saving and no29→28 theorem. A next useful bounded extension is to certify all second-root gates for these fixed witnesses, then identify whether a finite-duty classification forces a nonnative height gap uniformly over a genuinely complete pair of star domains.

## Verification

```sh
python3 scripts/fresh_compact29_coupled_first_roots_20261018.py --verify
python3 scripts/check_compact29_coupled_first_roots_saved_20261018.py
```

The producer uses the inherited exact coefficient field for projections and rational enclosures for labels11/12. The standalone checker imports no project modules: it reconstructs the integer E7/E8 Gram matrices, rational inverses, coefficient products, literal slot-pair quotient, strict signs, and independent tighter cosine enclosures. It does not repeat the original local geometric witness proofs or assert full-star/cone completeness. Generation f767fe, separate fresh replay60377b and standalone checker06ca8b all terminated with PASS. Metadata ee49bc passed7 source bindings,5 local links and2 ASTs.

Root completed a full mathematical producer/report/standalone-helper read: common finite-vertex sheet orientation, the exact11/12 native gap, and the complete49-root/819-quotient/3,087-literal declared domain are sound. The reviewer did not execute a program or reprove the inherited79 individual geometric witnesses. The packet is final/frozen with this precise independent review scope.
