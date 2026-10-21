# Endpoint and triangle domains at the expensive non-ADE vertices

The four non-ADE vertex types most prominent in the [small-D triangle census](fresh_sqrt5_rank13_small_d_third_types_20261017.md) now have complete necessary finite-endpoint and triangular-face domains. The [exact generator](../scripts/fresh_sqrt5_rank13_compensating_endpoints_triangles_20261017.py) and [saved output](fresh_sqrt5_rank13_compensating_endpoints_triangles_20261017.json) use rational arithmetic in Q(sqrt5).

The selected set is all rank13 forests with ADE, H3, H4 and I2(5) components whose inherited baseline is at least12. There are838 such forests in total and578 selected forests. Thus this extends the earlier215 selected ADE forests; it does not silently retain the ADE-only triangle exclusion.

| Original vertex forest | Baseline | Finite endpoint attachments | Triangle occurrences | All three corners selected | Possible original triangle pairs |
|---|---:|---:|---:|---:|---:|
| D4+D5+H4 |21|239|304|290|19|
| A4+D5+H4 |20|295|327|312|19|
| A1+2D4+H4 |19|372|495|474|26|
| 2D5+H3 |19|234|247|234|35|

The complete search tests23,036 endpoint contact assignments, retains1,140 finite endpoint attachments, tests6,840 framing labels and19,584 third corners, and retains1,373 marked triangle occurrences across142 third-vertex forests. Of these,1,310 have all three corners in the enlarged selected set. These are local necessary configurations, not counts of faces of any constructed polytope.

## Completeness and arithmetic gates

At a selected finite endpoint the new root has degree at most three. A trivalent node has only3-bonds; a node of degree two has at most one5-bond. The generator tests every subset and every such assignment, then identifies the whole destination forest by its complete weighted tree code. This includes isolated roots and changes of component decomposition.

All prescribed products are0,-1,-phi. If J is the inverse of the original finite Gram and i is the released node, put

    a=J_ii,
    b=-sum_j J_ij t_j,
    c=2-sum_jk t_j J_jk t_k,

where the new fixed contacts are -t_j. The remaining framing product g has Schur residual c-2bg-ag^2. The original Gram is positive definite at both real embeddings: after tree sign switches, conjugation reduces each absolute5-bond weight fromphi to1/phi, preserving positivity by comparison with the original positive-definite Stieltjes matrix.

The physical inverse is nonnegative, and its entries are strictly positive within each connected component. Therefore b is nonzero exactly when the prescribed contacts supply a nonzero path from the old port to the new root. Its cycle forces g into Q(sqrt5). Framing labels4,6,10 are then excluded. When b=0, their squares2,3,2+phi are sufficient to evaluate the residual. All six finite labels2,3,4,5,6,10 are tested, requiring a negative physical residual and a positive conjugate residual. No ADE rational-product exclusion is used.

For every retained edge and every other released original root, the resulting13-node third corner is tested against the full finite Coxeter list over the allowed alphabet, including the exceptional4/6/10 components. Its full forest, marked pair, baseline and exact Schur coefficients are saved. The accepted framing-label counts are2:663,3:368,4:22,5:279,6:4,10:37.

## Consequence for the next counting argument

The earlier rule that every triangle has a nonspecial vertex applied to an ADE selected set. It fails as an exclusion rule for this enlarged domain: the full local tests retain1,310 all-selected occurrences. Any extension of the [five-type incidence theorem](fresh_sqrt5_rank13_five_type_incidence_20261017.md) must allocate the required triangle weight using their actual three marked corner types. Assigning zero to these triangles would leave an unpaid obligation.

Generation and a fresh saved replay pass. A separate agent read the source and checked endpoint completeness, the path/cycle criterion, the conjugate-positive argument, Schur signs and third-corner bookkeeping. This was an internal mathematical/source review, not a separately implemented arithmetic verifier.

No dotted edge-length domain, quadrilateral census, simultaneous realization, new payment rule or dimension bound is claimed here.

```sh
python3 scripts/fresh_sqrt5_rank13_compensating_endpoints_triangles_20261017.py --verify
```
