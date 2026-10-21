# Arithmetic triangular faces at E7 and E8 over Q(√2)

The complete finite-contact one-root census leaves triangular completions at **10 of the 21 marked pairs in E7** and **9 of the 28 in E8**, for the arithmetic ground field Q(√2). The [exact checker](../scripts/fresh_arithmetic_sqrt2_local_20261015.py) and [data](fresh_arithmetic_sqrt2_local_20261015.json) use rational arithmetic and exact quadratic-field signs. These are local restrictions, not a dimension theorem; quadrilateral completions remain open.

Number the E_r diagram by the path `0—1—…—(r−2)` with the extra edge `2—(r−1)`. All old bonds have label3 and all root squares are2. The allowed marked pairs are:

| Component | Pairs admitting a triangular local star |
|---|---|
|E7|01, 02, 06, 12, 16, 23, 26, 34, 36, 45|
|E8|01, 02, 07, 12, 17, 27, 34, 45, 56|

For example, `06` denotes the unordered pair of nodes0 and6, not a bond label. All pairs absent from the table are excluded as triangular two-faces in the stated arithmetic setting. The argument allowing other orthogonal components of a full vertex diagram is given below. No assertion is made that the remaining local stars fit together in one compact polytope.

## Complete finite alphabet and squareclasses

Vinberg's arithmetic criterion places squared normalized pair products and all cyclic products in the ground field, with conjugate Gram matrices positive semidefinite. [Nikulin, §2, recalling Vinberg's criterion](https://arxiv.org/pdf/0708.3991).

For a finite angle π/m, `4cos²(π/m)=2+2cos(2π/m)` belongs to Q(√2). The real cyclotomic degree therefore implies `φ(m)≤4` for m≥3. The complete candidates are m=3,4,5,6,8,10,12. Direct evaluation removes5,10,12; including right angles leaves **2,3,4,6,8**. This is a field restriction, not a search cutoff.

Let p be the nonnegative contact row of the new root against the connected old E diagram. For any two nonzero contacts p_i,p_j, the old path joining i,j has products−1. Its cycle through the new root implies `p_i p_j∈Q(√2)`. Hence all nonzero contacts lie in one of three classes:

* mixtures of 1 and √2, corresponding to labels3 and4;
* only √3, corresponding to label6;
* only √(2+√2), corresponding to label8.

No mixed classes are missed. The first exceptional radical √3 is not in Q(√2); √6 is not either. The field norms of `2+√2`, `2(2+√2)` and `3(2+√2)` are2,8,18, respectively, so none is a square in Q(√2). These facts exclude every remaining cross-class product. Zero contacts are allowed everywhere.

The script enumerates all `3^r+2·2^r` profiles for each r=7,8, including duplicated zero rows. No dotted contact can occur in a triangular star: its new side shares a finite corner with each old black side, and every old white wall contains the face.

## Exact Schur tests and local realization

Let G be the old positive-definite Gram and J its rational inverse. The full one-root extension is Lorentzian exactly when

    q = pᵀJp > 2.

Its nonidentity conjugate is positive definite exactly when `q′<2`. Nondegeneracy makes strict positivity necessary even if one starts with the semidefinite version of the arithmetic criterion. For the label8 class, the two extensions of the field embedding differ by the sign of the new root and have the same Schur value.

Replacing old node i by the new root gives a finite vertex precisely when

    2 − q + (Jp)_i²/J_ii > 0.

Every pair of such replacement nodes is a possible marked pair of a triangular star. The original and the two replaced diagrams are positive definite; eliminating the white support gives a three-wall Lorentz Gram whose two-wall principal blocks are positive definite and off-diagonal entries nonpositive. Its three vertex rays are timelike in one cone and feasible, so they bound a compact triangle. As in the frozen [local-face construction](fresh_compact29_small_faces_20261015.md), a small inward white perturbation realizes the ambient local chamber. Conjugate positive definiteness rules out affine principal subdiagrams.

| Component | All profiles | Lorentz and conjugate-positive | Profiles with two or more finite replacements | Marked pairs covered |
|---|---:|---:|---:|---:|
|E7|2443|19|10|10|
|E8|7073|13|7|9|

The label6 class contributes no survivor. Label8 contributes two E7 profiles and one E8 profile, all with at least two finite replacements. The remaining survivors use labels3 and4. Cyclic products of these local stars are algebraic integers; the label8 star has only one nonzero contact and its squared product is `2+√2`. These local arithmetic checks do not establish finite volume of the ambient chamber.

## Other components of the vertex diagram

Suppose the full vertex diagram is E⊕Z, the marked pair lies in E, and x closes a triangular face. The subsets obtained by deleting either marked E node and adding x are spherical. If x is right to E, the full indefinite star would already be indefinite on Z+x, which is contained in both finite corners: impossible.

If the E+x Schur cost exceeds2, the census applies. Its two replacement tests are necessary because passing to a principal submatrix preserves positive definiteness. Its conjugate is positive definite because its determinant is nonzero and it is a principal submatrix of the conjugate ambient Gram. Extra components can only restrict the listed possibilities.

The exact inverse of E8 has minimum diagonal2 and strictly positive entries. Thus every nonzero contact row has cost at least2. Equality is an affine obstruction, excluded by compactness. For E7, the only nonzero profile with cost below2 is a single label3 contact at node5, with cost3/2. Multiple contacts already cost at least3 before their positive cross terms, and every other singleton has cost at least2.

In that E7 exception, x has at most one neighbor in any connected component of Z, since a second would produce a cycle in a finite corner. A corner retaining the E contact connects x to at least three nodes in every attached component. Such a finite connected Coxeter diagram cannot contain a label6 or8 bond. All its bonds therefore have label3 or4, whose squares are rational. Eliminating these finite trees contributes a rational Schur cost. Unattached components can be discarded. The total Schur cost is consequently rational and exceeds2, so its negative Schur complement persists at the other field embedding, contradicting arithmeticity. The affine cost-equality case is excluded as above.

This covers the potential loophole in which E+x alone is positive definite but x's contacts to other vertex components make the full star indefinite. It does not classify quadrilateral stars or combine different triangular stars globally.

## Status

Generation and saved replay pass. Independent review is recorded separately when complete. The script source and data are frozen; any stronger search uses a new namespace. The earlier unrestricted local-face witnesses remain valid in their stated nonarithmetic scope and are unchanged.

```sh
python scripts/fresh_arithmetic_sqrt2_local_20261015.py --verify
```
