# Nonzero A2 profiles for actual H9 H3 facets

There are exactly **seven possible source-contact profiles** for a norm-two
new facet of the actual H3 face on `S6={0,1,2,3,5,6}`: the all-right profile
and six nonzero profiles. The six were absent from the previous right-S6
searches. A bounded numerical pilot finds 2,553 candidates in this broader
domain, but its combined finite-bank incidence test finds no compact
closure. No new exact wall or finite-volume chamber is claimed.

## Exact, all-label facet-profile lemma

The retained support is `H4(0,1,2,3)+A2(5,6)`, verified from the independently
reviewed [old contact inventory](fresh_h9_persistent_cusp_screen_20261016.json).
All roots are positively normalized to square two. For any proposed actual
facet of this H3, its component orthogonal to S6 must have positive square.
Equivalently, its projection cost onto the positive-definite S6 is less
than two.

The H4 inverse Gram is entrywise positive and has diagonal

`14+6sqrt(5), 21+9sqrt(5), 10+4sqrt(5), 3+sqrt(5)`.

Every diagonal exceeds two. A nonzero acute Coxeter product has magnitude
at least one, so even a single nonright contact with H4 would already have
cost greater than two. All cross terms of the cost are nonnegative.
Consequently the new facet is right to the whole H4. This argument uses
exact inverse entries, not an imported finite-diagram classification.

The remaining inverse is

`G_A2^(-1) = (1/3)*[[2,1],[1,2]]`.

If both A2 contacts are nonzero, their magnitudes are at least one and the
projection cost is at least two. If a singleton contact has label at least
six, its magnitude is at least `sqrt(3)`, again giving cost at least two.
Parallel and dotted contacts are still larger. Thus only a singleton label
3, 4 or 5 on either node can occur, besides the zero profile.

| A2 contact | Positive residual square |
|---|---|
| Both right | `2` |
| One label3 | `4/3` |
| One label4 | `2/3` |
| One label5 | `(4-2phi)/3` |

Here `phi=(1+sqrt(5))/2`; the last value is positive. The
[source](../scripts/fresh_nonarithmetic_h3_a2_profiles_20261016.py) and
[certificate](fresh_nonarithmetic_h3_a2_profiles_20261016.json) check the
inverse matrices, positivity and exact threshold equalities symbolically.
This completes the possible profiles of an **actual facet of this H3**.
It does not classify arbitrary added ambient walls that miss the H3.

The compact target is necessary for a finite-volume ambient completion:
an ideal point of this face would be incident to its retained pi/5 pair.
The [persistent-cusp trace obstruction](fresh_h9_fixed_core_cusp_obstruction_20261013.md)
then forbids a full cusp translation lattice. This uses the entire retained
acute chamber and any common strict future interior, not an arithmeticity
assumption on new walls.

## Native offsets and the bounded endpoint scout

Write a new root as `r=p+y`, with `p` in A2 and `y` orthogonal to S6.
For prescribed A2 contact column `c`,

`p = (n5,n6)*G_A2^(-1)*c`,
`y^2 = 2-c^T*G_A2^(-1)*c`.

For every old wall `n_j`, the native product is

`(r,n_j) = (y,project_S6perp(n_j)) + (p,n_j)`.

The second term is retained for **all 55 old walls**. Likewise the common
center test uses `(p,H)+(y,H)<0`; it is not replaced by the zero-profile
time-coordinate inequality. For two new profiles, their native mutual
product includes `c^T*G_A2^(-1)*d` in addition to the residual Lorentz product.
Ignoring these offsets would produce false compatible caps.

The pilot uses the 13 old supporting edges with open ends, as reconstructed
in the [edge-continuation ledger](fresh_nonarithmetic_h3_edge_completion_20261016.md).
For each of six profiles and each edge, both endpoint contact labels range
from 2 through 10. The three projected normals must have positive-definite
Gram, so the new wall gives a finite endpoint. Of 6,318 recipes, 6,121 fail
this Schur test numerically. The remaining 394 conic branches yield 221
admissible open families, sampled separately at every original edge-end
coverage transition. Remaining nonconstant native products are strictly
dotted, with constant finite contacts allowed.

Additional isolated third-old-contact equations are solved separately,
without requiring their parameters to lie in the open all-other-dotted
intervals. There are 35,460 such equations. Together with the continuous
samples they produce 2,553 distinct numerical candidates. Some third-contact
solutions cut two or three original edge ends. Numerical single-cap clips
remain exterior; near-null output is unresolved and is not a certified cusp.

The domain retains the saved center. It does not cover arbitrary labels,
arbitrary first anchor pairs, arbitrary contact parameters or every possible
new-new coupled system. The exact seven-profile lemma must not be confused
with completeness of this numerical endpoint search.

## Combined finite-bank incidence test

The [joint source](../scripts/fresh_nonarithmetic_h3_a2_profiles_20261016_joint.py)
and [saved joint output](fresh_nonarithmetic_h3_a2_profiles_20261016_joint.json)
combine these candidates with the frozen zero-profile conic samples and
all 4,329 older three-contact recipes. The six trial walls used in some
older recipes are used only for reconstruction; they are not imposed as
additional retained walls. Native new-new offsets are included explicitly.

After deduplication and numerical norm checks, the pool contains 21,075
rows: 14,393 zero-profile samples, 2,353 nonzero-A2 candidates and the 4,329
old recipes. **3,023 sampled occurrences have uncertain norm residuals and
are skipped.** Several directions need large coordinates in the original
Lorentz chart, so this loss of precision is material to the search scope.

A finite intrinsic neighbor requires the residual two-by-two Gram to be
positive definite, in addition to native Coxeter compatibility. The pool
has 989 possible finite new-new edges. Every actual polygon facet of a
compact H3 must have at least three finite neighboring facets. Iterated
degree removal discards 19,469 and then 35 rows, leaving 1,571 candidates:
806 zero-profile, 485 nonzero-A2 and 280 from the old bank.

No surviving candidate cuts ideal-end target2 on old edge `[7,49]`.
Thus the degree-pruned finite model cannot close that target. The separate
279,814-row compatible-cover binary program is also infeasible. This is a
finite floating-point diagnostic, not a certificate against all continuous
profiles. It does not establish a fixed-wall noncompletion theorem.

## Next target and replay

The uncovered target lies in a very small Klein neighborhood where many
candidate normals have coordinates of size `10^4` or larger. The next
experiment should recenter the residual Lorentz chart near edge `[7,49]`,
reconstruct its native gates at higher precision, and solve coordinated
finite endpoint contacts there. The discarded norm-uncertain occurrences
must not be treated as geometric rejections.

Generation and fresh saved replay pass for the profile source; the joint
source is replayed independently as a saved numerical calculation. All
mathematical profile claims are exact; all conic, clipping, pair-graph and
binary-program conclusions are explicitly numerical. No frozen previous
source or result is modified.

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_nonarithmetic_h3_a2_profiles_20261016.py --verify
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_nonarithmetic_h3_a2_profiles_20261016_joint.py --verify
```
