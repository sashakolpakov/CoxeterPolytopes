# Simultaneous symmetric cusp filling in H9

This construction did not produce a compact H9 polytope. It instead gives
an exact obstruction for a continuous family that changes whole sets of
old right angles simultaneously. The obstruction applies to every Coxeter
label, not only the numerical alphabet used to discover it.

## The enlarged construction

The known 18-facet chamber `T2_D16_W1_5_9_14` has two center facets and eight
interchangeable pairs of leaf facets. In the
[source archive](../data/enumerations/alexeev_h9_doubling_cubes.json), the
centers are 1 and 3; the two leaf lists, paired entry by entry, are

```
u = [0,2,6,9,11,12,14,15],
v = [17,16,13,10,8,7,5,4].
```

Normalize the normals to square two. The old values are `a=b=0`, `w=2`,
and the first center has pairings −1 with the u leaves and zero with the
v leaves; the second center has the reverse pairings. The centers are
orthogonal.

The new family allows arbitrary nonnegative weights

```
(u_i,u_j) = (v_i,v_j) = −a       (i≠j),
(u_i,v_j) = −b                   (i≠j),
(u_i,v_i) = −w.
```

Any additional wall is allowed its own constant pairings −p with all u_i
and −q with all v_i, where p,q≥0. The extra walls need not have the same
p,q, nor be exchanged by a symmetry. All retained walls must admit a
strict timelike interior point, and their Gram matrix must have Lorentz
rank ten. Thus the new family changes 56 old right angles within the two
species and 56 more across different pairs. It is outside the preceding
[fixed-orthogonality obstruction](fresh_h9_cusp_filling_20261013.md).

**Theorem.** No compact Coxeter H9 chamber has this form. If at least one
of the two formerly right-angle weights a,b is nonzero, the declared
half-spaces have infinite volume. Arbitrary additional walls with the
constant pairings just described do not remove the obstruction. The
unchanged case `a=b=0,w=2` is excluded only from compactness, as expected
for the original finite-volume chamber.

This does not exclude nonsymmetric angle changes, deleting leaf facets,
or adding walls with nonconstant pairings inside a leaf species.

## Rank forces two branches

Take the fourteen differences
`u_i−u_8` and `v_i−v_8`, for i=1,…,7. Their Gram matrix is

```
[[2+a, b−w], [b−w, 2+a]]  tensor  (I7+J7).
```

The second factor is positive definite of rank seven. If the first factor
were nonsingular, this restriction would have rank fourteen,
contradicting ambient rank ten. Since 2+a>0, its determinant vanishes
precisely in the two branches

```
w = 2+a+b,       or       w = b−2−a.
```

The following explicit vectors work independently of the extra walls'
mutual pairings and without solving any remaining Gram determinants.

## First branch: a balanced leaf escape

Put `s=a+b`, `d=2+a`, and `S=Σ_i(u_i+v_i)`. For `s>0`, choose four signs
ε_i=+1 and four ε_i=−1, and define

```
Y = d S + 4s Σ_i ε_i(u_i−v_i).
```

Exact identities give

```
S² = −128s,
(S,Y) = −128sd,
Y² = 128sd(3a+4b−2),
(Y,u_i) = 8sd(−1+ε_i),
(Y,v_i) = 8sd(−1−ε_i),
(Y,extra wall) = −8d(p+q).
```

All wall pairings are nonpositive. A Coxeter weight for square-two
normals is either zero or at least one: it is `2cos(π/m)` when the walls
meet, two for parallel walls, or greater than two for disjoint walls.
Thus `s>0` implies `3a+4b−2>0`, so Y is spacelike. The timelike feasible
vector S has negative pairing with Y.

If s=0, this branch has w=2. Each matched pair is parallel and has a
positive affine kernel vector, preventing compactness under every
acute-compatible completion. No infinite-volume conclusion is made for
that unchanged branch.

## Second branch: an escape at the last pair

Here nonnegativity of w gives `b≥2+a`. Put

```
s=a+b,       d=2+a,       l=4−6a−8b,
Y = s S + l(u_8+v_8).
```

The identities are

```
S² = 16l,
(S,Y) = 4ld,
Y² = 4ld(l+s),
(Y,u_i) = (Y,v_i) = 0             (i≠8),
(Y,u_8) = (Y,v_8) = 2ld,
(Y,extra wall) = −2d(p+q).
```

Writing `h=b−2−a≥0` yields

```
l = −12−14a−8h < 0,
l+s = −10−12a−7h < 0.
```

Consequently Y is spacelike, is feasible for every wall, and has negative
pairing with the feasible timelike S. This branch has infinite volume
even without imposing a finite-angle alphabet on a or b.

For both branches, a strict timelike interior is assumed. If S lies on
some extra walls, perturb it slightly toward that strict interior point.
Then perturb Y slightly toward the resulting strict timelike point.
The resulting vector is strictly feasible and remains spacelike, with
the correct time orientation. A neighborhood of ideal directions lies
inside the chamber, proving infinite hyperbolic volume. The construction
does not rely on numerical vertex enumeration or on facet irredundancy.

## Numerical construction search

The preceding symbolic result was suggested by a numerical search with
two center walls exchanged together with u and v. Its finite-angle
alphabet is `m=2,3,4,5,6,8,10,12`, supplemented by seven declared dotted
weights. The rank-ten condition is imposed in either of the two remaining
two-by-two symmetric blocks. Ninety symmetry types cover all nine-facet
active subsets, including subsets of nonsimple vertices.

The scout examines 108,000 parameter tuples, with overlaps between rank
branches. It proposes 25,013 Gram matrices; 3,240 fail the numerical
signature gate and all 21,773 remaining proposals have numerical escapes.
No finite-volume candidate survives. These counts are numerical and not
counts of distinct chambers. Degenerate rank equations are not completely
enumerated by the scout. The exact theorem above covers them independently.

## Breaking the symmetry with seventy exact caps

The obstruction concerns extra walls fixed by the leaf permutations. To
leave that class constructively, specialize to the golden field
`phi²=phi+1`, with

```
a=0, b=(1+2phi)/2, w=2+b,
p=phi², q=0, c=2phi².
```

Here the two center walls have the exchanged pairings used in the numerical
scout, and their mutual pairing is −c. The 18-wall seed has signature
`(9,1,8 zero)`; its conjugate is positive semidefinite of rank ten.
All seed pairs are orthogonal or disjoint. This is a new seed in this
search, not a claimed new finite-volume polytope. Conjugate positivity
does not establish arithmeticity: the cyclic product
`b²=(5+8phi)/4` is not an algebraic integer.

For each four-element subset of the eight leaf pairs, give its indices
sign +1 and the others sign −1. The normal

```
r_epsilon = (2−phi) S/8 + phi Σ_i epsilon_i(u_i−v_i)/8
```

has square two. Its leaf pairings are zero and −phi, hence angles π/2 and
π/5, and its pairing with each center is −1, hence angle π/3. Every one
of these **70 normals is individually compatible with all 18 seed walls**
and cuts its corresponding balanced escape. These caps have nonconstant
leaf pairings, so the all-parameter obstruction above does not apply to
their completed configurations.

The interaction between caps is restrictive. If two sign patterns have
Hamming distance h, then

```
(r_epsilon,r_eta) = 2−phi² h/4.
```

For different balanced patterns, h is 2, 4, 6, or 8. Equal-square-two
normals with pairing in `Q(phi)` can have finite Coxeter pairing only
0, −1, or −phi. This follows from the degree of `2cos(pi/m)`: degree at
most two leaves m=2,3,4,5,6, and the square roots for m=4,6 are not in
`Q(phi)`. The displayed pair values pass the Coxeter gate only when h=8,
giving dotted pairing −2phi. Thus the 70-cap dictionary has exactly
35 compatible complementary pairs, and no compatible triples.

One complementary pair gives an exact 20-wall Coxeter cone. A balanced
escape with plus indices `{0,1,4,5}` survives the pair whose plus subsets
are `{0,1,2,3}` and `{4,5,6,7}`. Adding `S/1024` makes it strictly feasible
while it remains spacelike. The verifier transports and checks this
strict escape separately for all 35 complementary pairs. Every compatible
subset of this particular dictionary therefore still has infinite volume.
Other nonsymmetric mirrors on this new seed remain unrestricted.

A subsequent [complete leaf-angle census](fresh_h9_symmetric_finite_caps_20261013.md)
strengthens the dictionary result: every norm-two real normal whose
sixteen leaf pairings belong to `{0,−1,−phi}` and whose center relations
are Coxeter-compatible is one of these 70 caps. Arbitrary real center
pairings, including dotted distances, are allowed. The proof solves
30 assignment orbits exactly; other leaf labels and dotted leaf relations
remain outside its scope.

The [cap producer](../scripts/fresh_h9_symmetric_balanced_caps_20261013.py),
[exact archive](../data/enumerations/fresh_h9_symmetric_balanced_caps_20261013.json),
[independent Fraction verifier](../scripts/fresh_h9_symmetric_balanced_caps_verify_20261013.py),
and [PASS output](fresh_h9_symmetric_balanced_caps_verify_20261013.json)
record the Gram matrix, all coefficients, both signatures, all 2,415
cap-pair gates, and all 35 strict escapes. A forged cap is rejected.
An [independent agent review](fresh_h9_symmetric_balanced_caps_independent_review_20261013.md)
reran the verifier and checked the signature decomposition, finite-angle
alphabet, compatibility census, and all-subsets conclusion; no issue was found.

## Independent verification

The [symbolic producer](../scripts/fresh_h9_symmetric_filling_proof_20261013.py)
uses SymPy only to write polynomial identities. A separate
[standard-library verifier](../scripts/fresh_h9_symmetric_filling_verify_20261013.py)
reconstructs the generic Gram matrix over a rational polynomial ring,
checks every entry of the fourteen-contrast factorization, both matched
weight branches, both witness norms, all sixteen leaf pairings, a generic
extra-wall pairing, and the orientation identities. It rejects a corrupted
norm certificate. The generic extra-wall calculation applies separately
to every additional fixed wall. Result: **PASS**.
Another agent independently reran this verifier and reviewed the rank and
infinite-volume arguments, including the unchanged parallel-pair case.

- [Exact polynomial certificate](../data/enumerations/fresh_h9_symmetric_filling_proof_20261013.json)
- [Independent replay output](fresh_h9_symmetric_filling_verify_20261013.json)
- [Numerical scout](../scripts/fresh_h9_symmetric_filling_20261013.py)
- [Numerical archive](../data/enumerations/fresh_h9_symmetric_filling_20261013.json)

```sh
python scripts/fresh_h9_symmetric_filling_proof_20261013.py
python scripts/fresh_h9_symmetric_filling_verify_20261013.py
python scripts/fresh_h9_symmetric_filling_20261013.py
python scripts/fresh_h9_symmetric_balanced_caps_20261013.py
python scripts/fresh_h9_symmetric_balanced_caps_verify_20261013.py
```
