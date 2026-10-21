# Bounded white-contact continuation census

For each classical affine diagram of order 18, 19 or 23, the complete saved
one- and two-contact domain contains no white-contact configuration that both
continues by a new simple leaf and avoids all affine subdiagrams of order at
most 9. This is a bounded exact diagnostic, not an all-rank theorem.

The [source](../scripts/fresh_large_affine_white_continuation_scout_20261017.py)
and [certificate](fresh_large_affine_white_continuation_scout_20261017.json)
enumerate 21,492 contact profiles. They find 1,104 marked white-contact cases
with the two required spherical deletions. Of these, 411 retain both spherical
deletions after continuation. Every one of those 411 has an explicit small
affine witness already in the original enlarged diagram, before continuation.

| Affine kind | Order | White-contact marked cases | Continuing marked cases | Continuing cases without a small affine witness |
|---|---:|---:|---:|---:|
| A | 18 | 126 | 18 | 0 |
| B | 18 | 74 | 37 | 0 |
| C | 18 | 24 | 18 | 0 |
| D | 18 | 124 | 56 | 0 |
| A | 19 | 133 | 19 | 0 |
| B | 19 | 75 | 38 | 0 |
| C | 19 | 25 | 19 | 0 |
| D | 19 | 125 | 57 | 0 |
| A | 23 | 161 | 23 | 0 |
| B | 23 | 79 | 42 | 0 |
| C | 23 | 29 | 23 | 0 |
| D | 23 | 129 | 61 | 0 |

## Literal domain and checks

Start with the displayed labelled affine graph A of order N. Add z with
support of size one or two in A, with each nonzero contact having label 3,
4 or 5. For every pair of distinct marked roots b1,b2 of A, require both
`(A−bi)+z` to be spherical and require at least one contact of z outside
`{b1,b2}`. There are 417 contact profiles with at least one such marking;
the 1,104 count retains the different marked pairs.

Next add y, right to all of A and joined to z by label 3. Test both
`(A−bi)+z+y` for sphericity. The 411 continuing markings arise from 240
contact profiles. All of these continuing profiles have only label 3
contacts; the saved records retain their literal support and labels.

Sphericity is decided by exact labelled-tree isomorphism against the finite
Coxeter diagrams: A, B, D, E, F, H and the relevant dihedral cases. The
inherited [finite recognizer](../scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py)
contains every required rank and label in this bounded domain. Its tree
test rejects cycles before canonical coding. This uses the established
classification of finite Coxeter diagrams, not a numerical eigenvalue test.

For each continuing profile, the source searches connected induced subsets
containing z through order 9. A candidate must match an all-3 affine cycle
or a complete labelled affine tree. Its singularity and positivity are
then checked independently of that graph match: exact arithmetic over
Q(sqrt(2)) recovers a strictly positive kernel, checks a positive-definite
principal deletion, and verifies the connected weighted-Laplacian identity.
The literal roots, edges and kernel are saved. Thus every exclusion has a
direct positive-semidefinite corank-one certificate. No unreturned search
result is used as an exclusion.

The exact quadratic arithmetic comes from the separately frozen
[counterexample package](fresh_large_affine_white_contact_counterexamples_20261017.md).
All three source hashes are bound in the saved JSON. The two all-rank
counterexamples there remain valid: they show that contact only with marked
roots cannot be deduced from the two spherical deletions and the absence
of small affine diagrams alone. Their explicit extra-leaf continuation
failures are consistent with the census here.

## Scope and validation

The saved finite orders do not establish a conclusion for every N. Likewise,
the support-size and label ranges are the enumerated domain, not reductions
proved by this script. Applying the result to an arbitrary actual
quadrilateral requires separate all-rank arguments for the contact support,
labels, and continuation along the shared white path. No ambient chamber
realization, dimension bound, or completed quadrilateral-transfer theorem
is claimed here.

Generation and a fresh saved replay both passed, with zero unresolved
continuing markings. A separate narrow source review by the construction
agent checked the classical diagram constructors, deletion and white-contact
filters, finite recognizer use, and the exact affine-witness branch. That
review did not duplicate the census or arithmetic replay.

**FINAL/FROZEN:** source/data generation, fresh saved replay, source syntax,
all direct source hashes and local report links passed.
