# Independent review of the Q(√2) triangular-face restrictions

**PASS.** The [primary source](../scripts/fresh_arithmetic_sqrt2_local_20261015.py), [certificate](fresh_arithmetic_sqrt2_local_20261015.json) and [proof report](fresh_arithmetic_sqrt2_local_20261015.md) correctly give10 of21 possible triangular marked pairs in E₇ and9 of28 in E₈. The absent pairs remain excluded when the vertex diagram has additional orthogonal spherical components. This is an arithmetic triangular-face restriction, not a quadrilateral classification or dimension bound.

The [independent checker](../scripts/fresh_arithmetic_sqrt2_local_20261015_review.py) and [review certificate](fresh_arithmetic_sqrt2_local_20261015_review.json) import no producer or project arithmetic. They reconstruct the literal E diagrams from their three arms and invert their Cartan matrices with SymPy. Costs are evaluated by direct dot products in `QQ(sqrt(2))`; signs use shrinking rational intervals for√2 rather than the producer's squared-sign formula. Every accepted replacement is separately checked by a direct Gram determinant, without the producer's inverse-column replacement formula.

The independent enumeration reproduces all9,516 profiles, including the intentionally duplicated zero rows, all32 Lorentz/conjugate-positive survivors, their32 full Gram determinants, and237 replacement Gram determinants. It agrees with every recorded replacement sign and every covered or excluded marked pair. The complete subcritical list is also reconstructed:

* E₇: a single label3 at node5 has cost3/2; a single label3 at node0 has cost2.
* E₈: a single label3 at node6 has cost2; there is no nonzero profile with cost below2.

## Completeness of the contact domain

The proof's finite alphabet2,3,4,6,8 is sound. The real cyclotomic degree givesφ(m)≤4. Completeness of the small integer list can be checked without a cutoff assumption: φ(m)²/m is a product of prime-power factors `(p−1)² p^(e−2)`, all at least1 except the possible factor1/2 fromp=2,e=1. Hence m≤32, and the seven candidates are3,4,5,6,8,10,12. Their squared cosine fields remove5,10,12.

For two nonzero new contacts, their cycle through the rational old E tree puts their product inQ(√2). The three classes used in the producer are consequently complete. The norm obstructions2,8,18 exclude the products involving√(2+√2); elementary quadratic coefficient comparison excludes√3 and√6 fromQ(√2). Thus no mixed squareclass was omitted.

The Schur tests have the right strictness. When the physical one-root determinant is nonzero, its nonidentity conjugate cannot vanish, so conjugate semidefiniteness becomes positive definiteness. Deleting an old node leaves a positive-definite old principal block; the sign of the direct replacement determinant is therefore exactly the finite-corner test. The three finite corners and Lorentz three-dimensional residual define an actual compact triangular face of the local chamber. They do not prove compactness or finite volume of that chamber.

## Why other spherical components do not invalidate the exclusions

This was the main scope risk. Write the full vertex diagram as E⊕Z and let the new root x close a triangular two-face marked at distinct E nodes i,j. Both `(E−i)⊕Z+x` and `(E−j)⊕Z+x` are spherical. The full star is indefinite.

If x is right to E, all indefiniteness lies in Z+x, contradicting either finite corner. If its E cost is greater than2, the producer's E-only tests are necessary: both replacement matrices and the conjugate E+x matrix are principal submatrices of the corresponding full matrices. Additional components cannot repair a failed positive-definiteness test.

The exact subcritical census handles the remaining cases. A cost equal to2 is a nonzero positive-semidefinite one-root extension, hence an affine obstruction to compactness. The only case strictly below2 is the E₇ singleton3 contact at node5.

In that exception, choose the one of the two finite corners which retains node5. Every attached component of Z has at most one neighbor of x, because two neighbors would make a cycle in the finite diagram. Its connection to x and node5 makes a finite connected component of rank at least3. Such a component cannot contain a6 or8 edge: a high bond next to a nonright bond already has a nonpositive three-node path determinant. Consequently every edge of every attached Z tree and its contact to x has label3 or4.

Eliminating these trees gives rational Schur costs. This follows by leaf elimination, since squared edge entries are rational; equivalently, the relevant inverse diagonal and squared contact are rational. Unattached components, including any isolated high I₂ factors, are orthogonal and can be discarded. Adding these rational costs to3/2 makes the full cost rational. Since the full physical star is indefinite, that cost exceeds2. The same negative Schur complement persists under the nontrivial embedding ofQ(√2), contradicting arithmetic conjugate semidefiniteness. This proves that the subcritical E₇ profile cannot produce an omitted ambient triangular face.

Thus the exclusion direction is valid with arbitrary extra spherical components. The positive local stars remain only local possibilities; the review makes no assertion about simultaneous triangular stars or any quadrilateral.

## Replay and bindings

Generation and saved replay pass:

```sh
python scripts/fresh_arithmetic_sqrt2_local_20261015_review.py --verify
```

The review binds the primary PY, JSON **and MD**, as requested, along with its own checker. The primary report must therefore stay frozen; this review can be linked from a separate checkpoint. Published spherical-diagram classification, arithmetic cyclic-product/conjugate criteria and compact face correspondence remain theorem inputs reviewed at the proof level, not independently reproved here.
