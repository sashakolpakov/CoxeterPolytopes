# Same-type edges do not restrict the rank-21 vertex frequencies

**No dimension bound is improved.** The new result is an exact obstruction to a proposed refinement: arithmetic compatibility of the two endpoints of an edge, even with the deleted root specified, does not remove any of the finite vertex types in the current dimension-21 model. Every root port has a same-type V-arithmetic edge. The latest restricted dual above 40 therefore admits balanced endpoint/root-port flows.

The [checker](../scripts/fresh_sqrt5_universal_self_edges_20261016.py) and [certificate](fresh_sqrt5_universal_self_edges_20261016.json) use exact arithmetic in `Q(phi)`, where `phi²=phi+1`. They reuse the frozen finite component catalogue and the earlier exact field arithmetic, and check the recursive source manifests. This is local edge data and an incidence relaxation. It does **not** construct a simultaneous edge star, a common ambient lattice for all edges, valid two-face gluings, or a finite-volume polytope.

## A uniform construction

Let `A` be an integral symmetrized Gram matrix of a finite vertex diagram, positive definite at both embeddings. Choose root `i`, of norm `d=Aii`, and write

\[
s=\frac1{(A^{-1})_{ii}},\qquad
\lambda=\phi^8=13+21\phi,\qquad \tau=d\lambda.
\]

Duplicate root `i`, preserving all its contacts to the other roots. Give the two framing roots mutual product

\[
g=d(1-\lambda).
\]

Since `lambda>2`, their normalized product `g/d` is less than −1: the pair is dotted. Eliminating the common roots gives

\[
\begin{pmatrix}s&s-\tau\\s-\tau&s\end{pmatrix},
\]

whose eigenvalues are `tau` and `2s−tau`. The physical form has precisely one negative direction, because `s≤d` and `lambda>2`.

Conjugation sends `lambda` to `lambda⁻¹`. The conjugate block is positive definite whenever

\[
\lambda>\frac{d'(A^{-1})'_{ii}}2.
\]

The checker verifies all **724 root positions in 68 finite component types through rank 21**. The largest right-hand side is exactly **30**, attained in E8, whereas `phi⁸>30`. For this particular scaled rule, `phi⁶<30` fails the E8 test. Orthogonal spectator components preserve every check. A separate exact multiset count gives **311,574 rank-21 finite forests**, all covered by the same argument; this is not an extrapolation from sampled forests.

The finite Grams include the non-simply-laced types. Their root norms need not all equal two. For bonds 4, 6 and 10, integral symmetrizations use norm ratios `2`, `3` and `phi+2`, respectively. Every finite reflection coefficient `2Aij/Aii` is integral, as are the new framing coefficients `2(1−lambda)`. The dotted pair's length-two cyclic product

\[
4(1-\lambda)^2
\]

is integral and generates `Q(phi)`. Thus the construction supplies admissible first-type **V-arithmetic local edge data**, not a claim that the reflection subgroup has finite covolume.

## Actual endpoints and infinite-volume escape

These are actual finite hyperbolic edges, not just pairs of spherical diagrams. The checker computes the negative Gram-dual columns at the two framing walls. Their wall pairings are `−e_i` and `−e_new`; both have negative square and their mutual pairing is negative. They lie in the same time cone. Their positive segment has exactly the common walls active in its interior, and both endpoint diagrams equal the original finite diagram.

Writing `Delta=tau(2s−tau)<0`, the endpoint squares are `s/Delta`, their mutual pairing is `(tau−s)/Delta`, and

\[
\cosh(\text{edge length})=\tau/s-1>1.
\]

All these identities and the actual ray coefficients are checked for every rooted component. The complete rank-21 forest **D5 + 2E8** additionally has all **21** edge witnesses recorded with full endpoint rays and explicit strict null escapes.

There is also a general escape argument for every displayed rank-21 edge chamber. Delete any common root `j`. The remaining old-root diagram is finite; adjoining the duplicate has negative Schur residual because `tau>2d`. Consequently the relevant inverse-Gram diagonal is positive. A small positive perturbation of the endpoint midpoint gives a strict interior timelike vector `u`; a small positive `u`-perturbation of `−G⁻¹e_j` gives a strict interior spacelike vector `w`. The checker supplies exact perturbations for all 21 exceptional witnesses. Put

\[
U=u^2<0,\quad W=w^2>0,\quad V=(u,w),\quad
D=V^2-UW,\quad t=\frac{\sqrt D-V}{W}>0.
\]

Then `z=u+t w` is null, all its wall pairings are strictly negative, and `(z,u)<0`. An open patch of the ideal boundary lies inside the chamber, proving **infinite volume**. The radical is an explicit real scalar extension. No nonzero null vector over `Q(phi)` is asserted; that would contradict conjugate positive definiteness.

## The existing dual survives endpoint/root-port balance

The checker rechecks literal rational weak duality for the latest [A4 source-root model](fresh_sqrt5_dimension21_a4_source_orbits_20261016.md), whose exact lower bound is

\[
\frac{121786916150321}{3043000000000}>40.02199>40.
\]

Its 506 forest multipliers are nonnegative after changing the dual sign convention to vertex mass. At every one of the 21 labelled root ports of a forest with mass `mu`, assign self-edge mass `mu/2`. The same-port arithmetic edge above is allowed, and the endpoint balance is exactly `2(mu/2)=mu`. All **10,626** balance equations are checked and saved. Therefore this refinement alone cannot remove the forest frequencies supporting the existing obstruction. The old face-coefficient interpretation remains a dependency of the frozen producer/replay; this package independently checks the literal weak-duality arithmetic and the added port balances.

Self-port flow need not mean a loop in an abstract graph. A saved round-robin factorization of the connected simple graph K22 has 21 port colors, with exactly one edge of each color incident to every vertex. Every vertex may carry any one fixed rank-21 forest. This witnesses endpoint/port incidence only. It does not supply the triangles, face diagrams, global Gram matrix, or hyperbolic realization that a polytope would require.

## Where the global obstruction can enter

The local edge witnesses cannot be treated as independently composable. For each pair of root ports in D5 + 2E8, the checker tries the specific gluing that duplicates both roots and copies every cross contact. After eliminating the common rank-19 support, the symmetric block is

\[
\begin{pmatrix}2\alpha-\tau_i&2\beta\\2\beta&2\gamma-\tau_j\end{pmatrix},
\]

where `[[alpha,beta],[beta,gamma]]` is the finite two-root Schur complement. All **210** blocks are negative definite. Together with the positive antisymmetric directions, each proposed gluing has negative index two and cannot lie in a Lorentz space.

This rejects only that copied-contact gluing. It does **not** exclude other new-to-new products, shorter edge parameters, different endpoint identifications, or longer two-faces. The next meaningful computation is to leave the new-to-new product free and impose the full rank condition, Lorentz signature, conjugate positive semidefiniteness, and all four spherical corner conditions simultaneously. Endpoint-type matching alone contains none of that information.

Generation and saved replay passed. The source, certificate and report are frozen after those checks.

```sh
python3 scripts/fresh_sqrt5_universal_self_edges_20261016.py
python3 scripts/fresh_sqrt5_universal_self_edges_20261016.py --verify
```
