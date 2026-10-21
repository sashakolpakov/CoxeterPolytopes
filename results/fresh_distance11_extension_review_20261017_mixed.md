# Radical isolation of the mixed L22/L25 overlap

In each of the nine marked contexts below, the mixed L22/L25 core is orthogonal to every additional supporting white root whenever there is no marked white path of length at most 11. A preferred adjacent corner consequently has an isolated marked component of order at most **16**. This is a local lemma; the classification establishing that these are the surviving mixed contexts is outside this certificate, as is any new ambient dimension bound.

Use the literal diagrams from the [frozen quasi-Lannér census](fresh_high_rank_quasilanner_census_20261017.json). Number L22 by 0 through 8, and map the L25 nodes by

```
(0,1,2,3,4,5,6,7,8,9) -> (9,10,11,12,13,14,15,16,17,8).
```

Thus L22 node 8 and L25 node 9 are shared. The left black pair is one of \(\{1,2\},\{1,3\},\{2,3\}\); the right pair is the image of one of those three L25 pairs. These nine cases also cover the nine reverse orientations after exchanging the copies. Every additional cross bond between their exclusive parts creates a marked path of length at most 11, checked independently for all nine contexts. Hence, under the no-short-path hypothesis, there is no free cross contact in the following core Gram matrix.

## The actual linear relation

The [exact checker](../scripts/fresh_distance11_extension_review_20261017_mixed.py) deletes the two fork centers 0 and 9. The other 16 roots have a positive definite Gram matrix; the two-center Schur complement is

\[
 \begin{pmatrix}-3/20&-1/10\\-1/10&-1/15\end{pmatrix}.
\]

It is negative semidefinite of rank one. The 18-root Gram matrix therefore has inertia \((16,1,1)\). Its radical is spanned by

```
z = (-4,-2,-3,-2,-1,-3,-2,-1,0,6,3,4,2,5,4,3,2,1).
```

This Gram relation is an actual relation among the ambient normals. Indeed, the corresponding linear combination is orthogonal to the core span, which contains a timelike vector, and has square zero. The orthogonal complement of that timelike vector is positive definite, so the combination vanishes.

## Preferred corners and supporting contacts

Choose node 2 from \(\{1,2\}\), node 1 from \(\{1,3\}\), and node 2 from \(\{2,3\}\), on each copy. Delete the other marked root on each side. All nine preferred corners are positive semidefinite, and their marked roots lie in one connected component. There are three finite components of order 16, two finite components of order 15, and four affine components of order 15.

Every proposed extra supporting white belongs to the full adjacent-corner star. For each finite component of the preferred corner, positivity and the entrywise nonnegative inverse restrict a nonzero contact to inverse-diagonal entries at most 2, since every nonzero norm-two Coxeter contact has magnitude at least 1. For each affine component, its strictly positive kernel forces all extra contacts to that component to be zero. The omitted black roots are included separately in the allowed contact domain.

The certificate gives the exact inverse diagonals or positive kernel for every component and records the resulting domain. In every case:

* the zero-radical node 8 is excluded from that domain;
* each permitted negative-radical node is within distance at most 2 of the left black pair;
* each permitted positive-radical node is within distance at most 1 of the right black pair.

Let \(e\) be an additional supporting white, with nonpositive core contacts \(b_i=(e,r_i)\). The actual relation gives \(\sum_i z_i b_i=0\). If any contact is nonzero, the contacts must include both a negative-\(z\) node and a positive-\(z\) node. The resulting marked path through \(e\) has length at most

\[
 2+1+1+1=5,
\]

contrary to the hypothesis. Thus every additional supporting white is orthogonal to the entire core. There can be no indirect extension of the preferred marked component either, because the first connecting root would already have a nonzero core contact.

The [actual corner-factor statement](fresh_cusp_factor_payment_framework_20261017.md) supplies the same interpretation at ideal corners: the varying affine factor is contained in the supporting white diagram plus the two side roots. Hidden incident walls elsewhere in the cusp do not enlarge this marked component.

## Verification and scope

Generation and a fresh saved replay both passed. The [certificate](fresh_distance11_extension_review_20261017_mixed.json) records the exact radical and Schur matrix, all nine preferred-corner contact domains, positive definiteness of each supporting white diagram with one black root, and positive semidefiniteness of all 36 adjacent corners. It also checks every possible additional exclusive-part cross bond for a short path.

The checker uses the frozen standalone `Fraction` helpers from the [three-node-overlap certificate](fresh_distance11_extension_review_20261017.md), without numerical linear algebra or producer imports. Six direct source hashes, Python syntax and five report links passed. The exact-hub agent independently read the written proof and source and found no gap; this was a conceptual/source review without a duplicate arithmetic replay. This package makes no finite-volume realization or full distance-11 classification claim.
