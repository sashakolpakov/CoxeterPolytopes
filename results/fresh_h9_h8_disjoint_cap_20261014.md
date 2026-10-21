# A changed H9 seed with a compact H8 facet and three disjoint caps

The exact construction has **15 actual facets, exactly 44 finite vertices, and exactly 167 finite edges**. Three of those edges join the compact H8 facet to new finite vertices. A strict spacelike escape proves that this chamber still has infinite volume. This is a changed retained core: none of the twelve walls of the previous mixed-D4 chamber is retained.

The construction uses the literal [compact Bugaenko H8 source](../data/known_examples/bugaenko_compact_h8.json), whose 41 finite vertices are re-enumerated here. It does not claim a compact H9 example or a finite-volume lattice.

## Exact construction

Write `φ=(1+√5)/2`, `a=(√5−1)/4`, and `d=√5/4=a+1/4`. The H8 source has form `Q8=⟨−2φ⟩⊕E8` and eleven norm-two wall normals `u_i`. Use the H9 form

`Q = Q8 ⊕ ⟨d⟩`.

Retain the eleven perpendicular walls `U_i=(u_i,0)` and add the bottom wall `R=(0,…,0,1)`. Its facet is exactly the compact source. Let `v_0,v_32,v_33` be the three source vertex rays normalized to have minimum nonzero wall contact one. Each satisfies `v_i²=−a`; every nonzero contact `−(v_i,u_j)` is at least one. Add

`S_i=(v_i,−1)`, for `i=0,32,33`.

Their squares are `R²=d` and `S_i²=d−a=1/4`. Each new cap is right to the eight source walls at its selected vertex and strictly dotted to every other wall. In particular `(R,S_i)=−d`, and `d²>d/4`. The three caps are mutually strictly dotted.

The full 105-pair check gives 79 right pairs, eight pairs at `π/3`, two at `π/5`, and sixteen dotted pairs. There are no parallel pairs. The source center extended by final coordinate `−1` is strictly timelike and strict at all fifteen walls. Projecting it orthogonally to each wall gives an explicit timelike point strict at every other wall, proving all fifteen facets actual.

## Complete finite vertices and finite edges

The source replay examines all `C(11,8)=165` supports, including eight rank-deficient cases. A future timelike positive sum of all source normals bounds a projective slice, and the feasible support rays give exactly the archived 41 source vertices. Each has eight incident source walls and a positive-definite incident Gram matrix.

The changed chamber has these 41 bottom vertices `(V_j,0)` and the three top vertices

`P_i=(v_i,−a/d)`, with `P_i²=−a/(4d)<0`.

Every point is checked directly for future orientation, its exact nine-wall zero set, all remaining strict wall signs, and a positive-definite rank-nine incident Gram. All 44 projective points are distinct.

This list is complete. A finite vertex involving a cap can use only its eight finite-angle source neighbors, which gives precisely its displayed top vertex. A finite vertex involving `R` is one of the 41 source vertices. A support using only source roots cannot be positive definite of rank nine in the source Lorentz space. Dotted wall pairs cannot belong to a finite vertex.

The complete list of pairs of vertices sharing eight incident walls gives 167 finite edges: the 164 edges of the compact source facet and three vertical edges. In the certificate's vertex numbering, the new edges are `(0,41)`, `(32,42)`, and `(33,43)`. There remain `44·9−2·167=62` finite-corner edge germs without a second finite endpoint in this complete census. Their other ends and any components without finite vertices are not classified here; no bounded-link ideal-vertex claim is made.

## Remaining global escape

With the archived source center `H8`, the literal vector

`E=(v_0+v_32+v_33+H8/16,−10)`

has

`E²=−6837/256+6543√5/512>0`.

Its pairings with all fifteen wall normals are strictly negative, and its pairing with the strict future center is negative. Thus the chamber contains an open ideal-boundary region and has infinite volume. The three selected corner truncations leave a global gap; the finite-edge count alone does not establish a finite-volume completion.

## Host and duplication checks

The earlier [single-coordinate H8 lifts](fresh_h8_lift_h9_20261013.md) change an old endpoint wall and use extra coefficients `2φ` or `φ−1`; both retain a positive affine null-root combination. Here the eleven source directions are unchanged and the new coefficient is `d=√5/4`. The determinant ratios have field norms

`Norm(d/(2φ))=5/64`, and `Norm(d/(φ−1))=5/16`.

Neither is a rational square, so neither ratio is a square in `Q(√5)`. Since the forms have even dimension ten, this also distinguishes their determinant square classes under similarity. This checks the two concrete earlier lift hosts; it is not a classification of all quadratic hosts in the repository. A targeted search of the archived H8 lift, prism, and cap work found no previous copy of these three exceptional-ray caps with this coefficient.

Under `√5↦−√5`, the H8 form is positive definite while `d` becomes negative, so the changed host again has signature `(9,1)` at that embedding. More directly, the Cartan two-cycle of `R` and any `S_i` is

`4(R,S_i)²/(R²S_i²)=4√5`.

Its nonidentity conjugate is negative. Consequently any finite-covolume completion retaining these two walls fails the necessary conjugate positive-semidefiniteness condition for quasi-arithmeticity. This is a conditional obstruction for a future completion; the present infinite-volume chamber is not a lattice.

## Artifacts and verification

- [Producer](../scripts/fresh_h9_h8_disjoint_cap_20261014.py)
- [Exact certificate](fresh_h9_h8_disjoint_cap_20261014.json)
- [Frozen scalar arithmetic dependency](../scripts/fresh_h9_mixed_d4_coupling_20261014_review.py)

All geometry is exact in `Q(√5)`. Scalars use the inherited serialization `[[a,b],[0,0]]` for `a+b√5`; the extra field coefficients are zero. The certificate stores the form, all normals and pairings, actual facet points, all finite vertices and edges, the escape, and source hashes. The saved replay command is:

```sh
python scripts/fresh_h9_h8_disjoint_cap_20261014.py --verify
```

Generation, the saved replay, syntax checks, source hashes, and local report links all pass. The [independent review](fresh_h9_h8_disjoint_cap_20261014_review.md) reconstructs the geometry with a separate SymPy algebraic-field implementation, including the complete source support census and all 44 finite vertices.
