# Independent review of the distance-eight high disjoint-core lemma

**PASS.** The mathematical/source review found no gap in the stated local lemma. A separate exact certificate checker reconstructs the whole finite domain and verifies every saved witness using its own Fraction and Q(sqrt2) arithmetic, without importing the producer or its arithmetic helpers.

The [primary source](../scripts/fresh_distance8_high_disjoint_20261017.py), [certificate](fresh_distance8_high_disjoint_20261017.json) and [primary report](fresh_distance8_high_disjoint_20261017.md) exclude two disjoint minimal indefinite cores, both of orders6 through10, in an actual quadrilateral star with no opposite-black-pair white path of length at most8. This is not an overlap theorem, a supporting-extension theorem, or an ambient dimension bound.

The independent [checker](../scripts/fresh_distance8_high_disjoint_20261017_review.py) and [review certificate](fresh_distance8_high_disjoint_20261017_review.json) inherit the listed high-core catalogue as their input. They do not reproduce the proof that this catalogue contains every possible high core. They independently verify the actual normalized Gram inertia of all26 listed cores, all142 listed finite deletions and all49 listed ideal deletions, including graphs whose label4 cycles prevent a rational positive diagonal gauge.

## Coverage and label reduction

The checker rebuilds all384 marked contexts and their graph distances, with maximum marked radius7. A cross contact touching a black root would produce a path of length at most8; every admissible cross contact therefore joins white ports. A contact between ports of marked distances d,e requires d+e≥8.

For each port, both black-deletion Gram matrices are positive definite. Their normalized inverse diagonals are reconstructed independently. Adding a new norm-two root with a nonright contact of magnitude beta gives the necessary strict bound `beta² Jii<2`. Since a finite positive-definite Z-matrix has nonnegative inverse, other acute contacts cannot weaken this bound. This leaves754 of2,404 ports.

For the two supporting-white blocks the necessary bound is `beta² wi wj<1`. It is again monotone under additional acute contacts, as can be seen from the nonnegative vectors obtained by applying each positive inverse to the contacted coordinate. Of12,135 eligible single joins,3,037 fail white positive definiteness already for a simple bond.

On every remaining join the independently checked rational bounds either exclude all labels at least4, or exclude all labels at least5. The latter uses `beta²≥(3+sqrt5)/2>5/2`; it also excludes parallel and dotted contacts. Exactly178 joins retain label4 in addition to label3. Thus these are exhaustive Coxeter-label reductions, not imposed finite-label cutoffs.

The zero-cross-edge case is immediate: the two disjoint indefinite principal blocks already supply two negative directions. The source enumerates every nonempty cross-edge subset surviving the necessary single-edge filters. There are10,990 subsets with at least two edges. The independent checker reconstructs this exact set and all remaining label products, rather than only verifying the number of saved records.

In particular, the eligible edge lists are populated **before** the full single-join inertia and corner rejections. A two-negative single-join matrix does not by itself justify removing that edge from configurations with additional cross contacts. The source correctly retains it and checks the multiple-contact configurations separately.

## Single joins and the zero-cofactor argument

The full high cores are nonsingular of negative index one. If deletion of a contacted node is affine, its principal cofactor is zero, so the corresponding full inverse diagonal is exactly zero. Eliminating that core across the single cross bond subtracts

`beta² (G_core inverse)ii e_j e_j^T = 0`

from the other core. The joined matrix has two negative directions. This is a real symmetric Schur-complement identity; it needs no rational Gram gauge for the full core. The independent checker confirms the affine cofactor, full-core inertia and finite deletions separately in exact Q(sqrt2).

There are262 saved zero-cofactor type certificates and four types with explicit rational inertia calculations. The type keys are unordered endpoint pairs; a stored explicit Gram may retain either original block order. The checker accepts only the corresponding literal block swap and verifies the actual gauge entries before calculating inertia. This does not identify different contact positions.

The resulting9,256 single-label cases with two negative directions are impossible in a Lorentz space. The other20 single-label cases are excluded by explicit **strictly negative** vectors in an actual corner principal matrix. Strictness matters: ideal corners may be positive semidefinite. The producer's conversion from a nonnegative zero-square vector v with `(Gv)_k<0` to `v−(Gv)_k/Gkk e_k` preserves nonnegativity and makes the square strictly negative. The reviewer directly checks each final saved vector and its negative square.

## Multiple contacts

For10,987 multiple-edge subsets, a nonnegative nonpositive-square vector already contradicts actual-edge positive definiteness when all selected cross bonds are simple. Keeping the same positive gauge and increasing any cross magnitude only decreases this quadratic form. The unit-bond fallback, where used, is also a valid relaxation: all actual nonzero acute products have magnitude at least one. The independent checker verifies the literal native matrix, gauge, vector coordinates and nonpositive square for every such record.

The three remaining subsets use two copies of L25 and cross positions `(3,9)` and `(9,3)`:

| Black pair on each copy | Remaining label assignments | Exact rejection |
|---|---|---|
| `{0,1}` | `(3,3)` | Inertia `(18,2,0)` |
| `{1,2}` | `(3,3)` | Inertia `(18,2,0)` |
| `{1,2}` | `(3,4)`, `(4,3)`, `(4,4)` | Actual-edge negative vectors, squares `−5/6`, `−104/77`, `−12/7` |
| `{1,4}` | `(3,3)` | Inertia `(18,2,0)` |

The omitted doubled assignments in the first and last rows have already failed the individual label gates. They are not missing cases. All six surviving assignments are independently checked.

## Verification scope

Independent generation and a fresh saved replay pass. The certificate records source/input hashes and all check counts. The exact checker uses no SymPy, numerical linear algebra, producer imports, or inherited arithmetic implementations. Its inputs are the high-core catalogue and the primary saved certificate; the geometric facts that actual edges have positive-definite supporting Grams, corners have positive-semidefinite Grams, and roots lie in one Lorentz space are mathematical premises of the local statement.

The review confirms complete high-core **disjoint** coverage and its exact obstructions. It neither promotes a new dimension theorem nor closes the separate overlap and supporting-extension obligations.
