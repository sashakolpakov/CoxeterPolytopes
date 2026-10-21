# Independent audit of the 728 endpoint caps

**PASS.** All 728 archived entries are distinct norm-two mirrors, with the prescribed Coxeter or strictly dotted relations to the 16 seed walls and a common strict timelike center. Exactly 32 mirrors occur in the earlier 3,040-entry bank; 696 do not occur in that literal bank. This is a comparison of two saved banks, not a claim of global novelty or a new finite-volume Coxeter polytope.

The [independent verifier](../scripts/fresh_h9_split_endpoint_families_independent_20261013.py) imports no producer. It uses the earlier independent Fraction implementation of Q(phi), exact signs in quadratic extensions, and integer-square-root enclosures whose endpoints are independently checked by exact field signs. Its [saved result](fresh_h9_split_endpoint_families_independent_20261013.json) records source hashes, exact counts, all overlap IDs, and all ten caps cutting the old escape. The audited input is the [728-cap exact archive](../data/enumerations/fresh_h9_split_endpoint_families_20261013_exact.json).

| Check | Exact result |
|---|---:|
| Normal-to-covector identities and norm two | 728 |
| Strict common-center inequalities | 728 |
| Spherical extension of an orthogonal rank-eight seed block to rank nine | 728 |
| Right-angle seed relations | 4,336 |
| pi/3 seed relations | 1,328 |
| pi/5 seed relations | 1,232 |
| Strictly dotted seed relations | 4,752 |
| Parallel seed relations | 0 |
| Single-cap affine subdiagrams with the seed | 0 |
| Distinct mirrors | 728 |
| Mirrors also present in the old bank | 32 |
| Mirrors outside the old bank | 696 |

The seed Gram matrix is rebuilt from a=phi+1/2, b=−a/(12a+1)=−5/22+phi/11. An explicit congruence basis verifies inertia (9,1,6); quotienting the coefficient kernel gives Lorentzian dimension ten. The all-ones coefficient vector has square −64−128phi and pairs strictly negatively with every seed and cap.

For every cap, the verifier multiplies its archived coefficients by this independently rebuilt Gram matrix, computes its square, and checks every seed pairing exactly. It also compares the resulting covector with rational intervals of radius 10^−6 around the original decimal proposal, independently selecting the correct real radical branch. For the chosen orthogonal eight-wall seed block, the Schur complement

```
q = 2 − (1/2) sum_i (cap, seed_i)^2
```

equals the archived value and is strictly positive. Thus the nine-wall extension is spherical. This alone does not assert that all these walls form a vertex of a larger chosen cap configuration; other retained inequalities must also be checked there.

The single-cap affine check is elementary. Every cross pair between the two seed blocks is strictly dotted, so a positive-semidefinite subdiagram cannot meet both blocks. Within one block, a cap and a subset of seed walls form a star. Strictly dotted arms cannot occur in a positive-semidefinite subdiagram; a finite-label star would be singular only if k3+k5 phi²=4. Irrationality of phi forces k5=0 and k3=4. Each archived cap has fewer than four pi/3 neighbors in either block, and has no parallel seed pair. This excludes affine subdiagrams using one cap and any seed subset. It says nothing about affine subdiagrams using two or more caps.

The old universal escape has coefficients

```
[8,8,8,1,-3,-3,-3,-3,-3,-3,-3,5,5,5,5,5].
```

Its square is (7242−2800phi)/11>0 and its center pairing is −116−232phi<0. It is strict on the seed. Its cap pairings are positive at indices 11, 21, 31, 92, 98, 108, 114, 202, 208 and 521, and negative at all other 718 caps. Indices start at zero. Hence this particular old escape does not exclude every selection from the new bank; no finite-volume conclusion follows from this sign test.

Distinctness and overlap use complete pair comparisons, not a numerical nearest-neighbor query. Certified covector intervals at scale 10^24 separate every unequal pair within the new bank and between the new and old banks. The 32 remaining cross-bank comparisons are checked as exact equalities in Q(phi,sqrt(D),sqrt(E)), allowing dependent radicals. Norm two and the common negative center orientation mean covector equality is exactly equality of oriented mirrors. The replay covers 2,477,748 pairs and checks 337 distinct coefficient enclosures; it finds no duplicate within the new bank. Every matched old covector is also recomputed from its saved normal. The old source is hash-anchored to the previously audited bank.

Reproduction from the repository root:

```sh
python scripts/fresh_h9_split_endpoint_families_independent_20261013.py --verify
```

The verifier rejects altered normal coefficients, a forged finite seed label and an incorrect escape-pairing sign. Both artifact construction and replay pass. The audit establishes individual caps, seed geometry and literal bank overlap; it does not certify cap-cap compatibility, completeness of the endpoint search, compactness, or finite volume.
