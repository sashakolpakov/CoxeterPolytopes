# Intrinsic endpoint profiles for the 22-wall cusp seed

The [exact census](fresh_h20_intrinsic_caps_20261014_endpoint.json) supplies every necessary finite or ideal Coxeter attachment to the nineteen walls of each of the seed's **252 outgoing cusp edges**. It allows arbitrary real new normals, including labels 5 and the H3/H4 diagrams. There is no ambient integrality or rationality restriction.

The literal seed is **affine A11 × affine A2 × affine E6**. The rank-six factor was previously called D6 in some descriptive messages, but its seven-node Gram is the affine E6 tree: center22 joins33,40,45, which continue to26,34,24. Its marks are `3,2,2,2,1,1,1`. The earlier positive-kernel and rank arguments do not depend on that mistaken name.

After normalizing each wall to square two, the three finite edge types and their complete contact lists are:

| Deleted E6 node | Old spherical support | Literal edges | Finite profiles | Semidefinite profiles |
|---|---|---:|---:|---:|
| Outer leaf | A11 + A2 + E6 | 108 | 19 | 9 |
| Middle node | A11 + A2 + A5 + A1 | 108 | 54 | 14 |
| Center | A11 + four A2 components | 36 | 71 | 93 |

There is also one isolated-A1 contact family in the middle type, containing every integer label `m≥7`, and its parallel limit. Labels 3–6 on that wall already occur in the finite lists. The three finite lists contain 260 profiles, including the all-right profile in each type. Transport to all literal edges gives 10,440 finite and 5,832 semidefinite profile occurrences, plus 108 arbitrary-label families. These are local attachments, not a mutually compatible cap selection or a finite-volume construction.

## Complete endpoint test without an arithmetic angle alphabet

Let B be the square-two positive-definite Gram of the nineteen retained seed walls. Write the prescribed cap pairings as `−c_i`, where each c_i is either zero or `2cos(π/m)`, and normalize the cap to square two. The full twenty-wall Gram is positive semidefinite precisely when

`cᵀ B⁻¹ c ≤ 2`.

Strict inequality gives a spherical endpoint diagram. Equality is retained: a proposed ideal endpoint cannot be discarded merely because it is not positive definite.

The old Gram is a positive-definite matrix with nonpositive off-diagonal entries, so its inverse is entrywise nonnegative; every inverse diagonal is at least `1/2`. The checker verifies these facts for all three literal representatives. Each nonright finite contact has c_i≥1, hence the Schur cost is at least half the contact-support size. There are at most four nonright contacts.

Suppose a contact has `m≥7`. If its old wall has an old neighbor, the three-node principal Gram contains a contact of magnitude greater than `√3` and an old contact of magnitude one. Its determinant is negative, even before any further nonpositive pairing is included. The same conclusion holds if the cap has a second nonright contact. Consequently such a label is possible only as the sole contact with an isolated old A1 node. A parallel contact also requires precisely that isolated two-node block. These are the explicitly saved tails.

Every other contact therefore has label 3, 4, 5, or 6. The producer enumerates all supports of sizes zero through four, all assignments from that pool, and the exact Schur inequality. The minimal-label support filter is safe because every inverse entry is nonnegative. All computations lie in the multiquadratic field with basis `1,√2,√3,√6,√5,√10,√15,√30`; exact zero detection and rational square-root enclosures decide the signs. No finite-type classification table or original-source angle theorem is needed for completeness.

The 252 transports are explicit permutations. Rotations of the A11 cycle, cyclic permutations of the A2 triangle, and permutations of the three E6 branches preserve every normalized Gram entry and positive affine mark. The certificate checks those identities rather than assuming orbit equivalence from the names of the diagrams.

## Compatibility with all 22 seed walls

Only the three omitted walls remain to test after fixing an endpoint profile. In the upper-half-space model take the product of the three Euclidean simplexes with inradii `1,s,1`, where s is positive. The original record0 at `u=1/2` has `s≈1.93786809257`. Exactly, `s²=8qL²/(k23² n23)`, where `n23=q*r23²−a23²` in the frozen quartic notation. Allowing other positive shifts changes s proportionally; the endpoint diagram census itself is unchanged. Let z be the product vertex for the chosen cusp edge, let n_i be the nineteen unit outward normals there, and solve

`n_i·y = cos(π/m_i)`, with zero for a right contact.

The saved normalized spatial cost is `y² = cᵀB⁻¹c/2`. A sphere cap of radius R has center

`w = z − R y`.

For each omitted wall o, let m be its positive affine marks. Its signed inward distance from w divided by R is exactly

`A_o/R − B_o`, where `A_o = r_factor (Σm_i)/m_o` and `B_o = Σ_retained (m_i/m_o) cos(π/m_i)`.

Thus `A_o` is 12 for A11, `3s` for A2, and `12/m_o` for E6. The checker saves all three B_o values in the exact number basis, along with numerical copies for proposing constructions. Compatibility with all 22 seed walls is equivalent to requiring each of these three ratios to be a finite Coxeter cosine or at least one, the parallel/divergent range. A negative ratio is obtuse and is rejected.

This formula has no hidden continuous variables besides R and the chosen positive triangle inradius s. In particular every finite profile has an explicit nonempty all-divergent range: choose R positive and strictly smaller than every `A_o/(B_o+1)`. The cap then retains all seed walls and cuts the chosen vertical edge at height `R√(1−y²)>0`. The twenty incident normals are independent and positive definite, and the other three seed inequalities are strict there, so this is an actual finite endpoint. This establishes intrinsic cap existence for each listed finite profile; it does not establish their mutual compatibility or global coverage.

For the semidefinite profiles, the point is at height zero. Their Gram test alone does not prove a bounded cusp link: additional incident walls may be required. None is recorded as a completed finite-volume ideal endpoint.

Equivalently, in Lorentz coordinates let D be the cusp null ray, U the nineteen-wall positive span, and V the future null vector in U-perp with `D·V=−1`. With `h=−r·D>0`, a cap with the prescribed pairings is

`r = −Σ(B⁻¹c)_i s_i + hV − (2−cᵀB⁻¹c)D/(2h)`.

This identity proves square two and every prescribed endpoint pairing. The three omitted-wall pairings are `(Σ_retained m_i c_i − μh)/m_o` when the normalized affine relation is `Σm_i s_i=μD`.

## Noncrystallographic constructive profiles

The outer type has two label-five profiles: attaching to either endpoint of its A2 component creates H3. The central type has eight such H3 profiles, one for each endpoint of its four A2 components. The middle type has two H3 profiles, one isolated I2(5) profile, and two H4 profiles.

The H4 profiles are middle-type entries43 and51. They contact isolated wall26 at label5 and respectively wall50 or200 at label3, giving the chain with labels `5–3–3`. Their normalized spatial cost is approximately `0.9878418305<1`. They give genuine finite endpoints for sufficiently small R, and are available to the separate joint-cap construction search. The original source's crystallographic-angle restriction does not apply after introducing these intrinsic walls.

## Data interface and verification

In each JSON `types` entry, `edge_support` specifies the order of the nineteen-entry `cosine_vector` and `cosine_vector_approx`. Each profile includes its sparse contacts, exact Schur cost, finite/semidefinite status, and three omitted-wall ratio functions. `literal_edge_transports` maps these vectors to every literal cusp edge without changing their order.

```sh
python scripts/fresh_h20_intrinsic_caps_20261014_endpoint.py --verify
```

Generation and saved-result replay passed. The [producer](../scripts/fresh_h20_intrinsic_caps_20261014_endpoint.py) binds the frozen seed and raw source. The current positive result is a complete supply of locally and seed-compatible endpoint constructions, including H3/H4; no assertion of a bounded H20 chamber is made.

The [separate review](fresh_h20_intrinsic_caps_20261014_review.md) passed, independently recomputing all 260 accepted profiles, all omitted-wall formulas, and all 252 transports. It source-reviews the exhaustive rejected-profile argument rather than claiming a second enumeration. The subsequent [all-radius bound](fresh_h20_intrinsic_caps_20261014_endpoint_bounds.md) shows why caps meeting other faces are necessary at the golden triangle scale.
