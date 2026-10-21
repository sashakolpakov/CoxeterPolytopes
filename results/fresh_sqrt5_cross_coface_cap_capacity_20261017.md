# A cap can occur only once across these orthogonal cofaces

The [exact certificate](fresh_sqrt5_cross_coface_cap_capacity_20261017.json)
and [checker](../scripts/fresh_sqrt5_cross_coface_cap_capacity_20261017.py)
give 241 additional signed capacities in the saved rank-13 arithmetic
model. Seventy have positive load at the current restricted face masses,
and 53 are violated. These constraints couple different immediate
cofaces of the same maximal simplex. Thus the old ratio
12.344944732908832 is not a dual witness for this enlarged model.
This certificate does not compute a repaired witness, a new optimum or
a dimension bound.

## Geometric capacity

Normalize every original outward wall normal to squared norm 2. The full
star of a maximal simplex S has a nonsingular Lorentzian core and an
orthogonal finite spectator forest. Consider only the exported prism-base
occurrences satisfying all of the following restrictions:

* The immediate coface root q lies in an isolated A1 spectator component.
* The complete coloured core automorphism group is trivial.
* Deleting the prism's other base normal b from its raw connected star
  leaves precisely the core and q. Every further spectator is orthogonal
  to b.
* The squared normalized b,q product is tau with 4 < tau < 6.

Fix the complete normalized contacts of b with the canonically labelled
core, its base-vertex subset B, and tau. Since the core is nonsingular,
these contacts determine a unique projection v into its span. The full
simplex star spans the ambient Lorentz space. Orthogonality to all other
spectators and the outward contact with q therefore force

    b = v - sqrt(tau) q / 2,
    (v,v) = 2 - tau/2,     -1 < (v,v) < 0.

The same profile placed at a different isolated spectator q' has normal
b' = v - sqrt(tau) q'/2. Since q and q' are orthogonal and both orthogonal
to the core, (b,b') = (v,v). This is impossible for two distinct walls of
a hyperbolic Coxeter polytope: a normalized nonzero negative wall product
has magnitude at least 1. The alternatives are 0, -2 cos(pi/m) for
m >= 3, -2, and products below -2.

At the same q, the normal is unique. With B fixed it contributes the same
base occurrence, not a second one; the inherited actual vertex/coface
flag capacity gives the same conclusion. Consequently, among all
isolated A1 copies of this S-star, the total number of selected occurrences
with this profile is at most one. The trivial core automorphism restriction
prevents any unrecorded change of the contacts under a core symmetry.

The resulting signed counting column is

    sum(selected actual prism-base incidences) <= number of S-faces.

Omitted types have zero positive coefficient. This weakens the valid
inequality and does not require classifying every possible cap. Different
columns are separate capacities; their right sides are not added as
disjoint resources.

## Exact reconstruction and the nine-copy violation

The checker reconstructs all eligible entries from the frozen raw stars,
actual provenances, base maps and labelled coface data. Of 1,827 isolated
A1 base occurrences, 847 have trivial core symmetry, 434 also have the
required sole internal A1 extra, and 241 have the exact dotted parameter
4 < tau < 6. They form 241 groups in this data set. All field comparisons
use rational arithmetic in Q(sqrt(5)); no angle tolerance is used.

For the particularly simple group with base type 54999, the S-star has
nine isolated A1 copies and a five-node Lorentzian core. Its selected prism
is actual provenance 440, merged face 26259, raw star 23. Its other base
has a different type and is not counted here. The relevant exact values are

    tau = (7 + sqrt(5))/2,
    (v,v) = (1 - sqrt(5))/4,
    alpha[26259] = 7551/10000000000,
    alpha[80143] =  839/10000000000.

The old face masses saturate all nine separate coface capacities, but
the new single capacity has load nine times its right side. The
[independent witness check](review_sqrt5_cross_coface_full_base_witness_20261017.md)
reconstructs the normalized core Gram matrix, its inverse and determinant,
the full base mapping and all core permutations. In particular it verifies
the projection norm directly, independently of the above norm identity.

## Verification and limits

Generation and fresh saved-data replay pass all exact reconstruction and
load checks. The
certificate binds its source and frozen inputs and invokes no optimizer.
The existing complete within-one-coface packing cones remain valid, but
they do not impose this cross-coface wall compatibility.

An enlarged-model witness must be repaired or recomputed before it can
be used. Raising simplex masses also changes vertex and local packing
charges; the old covers cannot be reused unchanged. Neither the failure
of the old witness nor a conservative repaired ratio below 12 establishes
primal feasibility. The complete unrestricted dimension-163 proof and
the compact arithmetic dimension-21 proof are unaffected.
