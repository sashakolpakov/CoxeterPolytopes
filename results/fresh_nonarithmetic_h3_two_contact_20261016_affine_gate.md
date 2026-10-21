# Exact exclusion of the two-contact null-line recipes

The [checker](../scripts/fresh_nonarithmetic_h3_two_contact_20261016_affine_gate.py)
and [certificate](fresh_nonarithmetic_h3_two_contact_20261016_affine_gate.json)
exclude **all real parameters in23 contact recipes**, accounting for all25
null-line branches in the new H3 pilot. Generation and saved replay pass.
This is an exact family obstruction, separate from the numerical searches.

Retain the old walls and let a proposed norm-two wall be right to
`S6={0,1,2,3,5,6}`. The retained pair0,1 has angle pi/5. Each of the following
anchor pairs is right and both anchors are right to0,1:

```
(8,44), (8,52), (30,33), (30,44),
(30,47), (34,47), (34,52).
```

For each pair, new contacts `(3,6)`, `(4,4)`, or `(6,3)` form an affine
triangle. In addition, the old pair `(30,52)` has angle pi/4. New contacts
`(2,4)` or `(4,2)` again form an affine triangle. These are21+2 distinct
recipes; the last two have two numerical null-line orientations each.

All old labels and right products are inherited from the independently
reviewed [49-wall contact inventory](fresh_h9_persistent_cusp_screen_20261016.json)
and [wall52 certificate](fresh_nonarithmetic_escape52_20261016_cap.json), with
its [independent review](fresh_nonarithmetic_escape52_20261016_review.md).
No numerical near-zero test is used to establish these old incidences.

The checker reconstructs every three-root norm-two Gram using the exact
values `0,1,sqrt(2),sqrt(3)` for `2cos(pi/m)`. It verifies rank2, all positive
two-by-two principal minors, and a strictly positive null kernel. For example,
when the old anchors are right, the kernel is proportional to
`(2cos(pi/m1),2cos(pi/m2),2)`.

Let nu be the resulting positive linear combination of the three roots.
In any coherent acute realization sharing a strict future interior point X,
one has `nu.X<0`. Thus nu is nonzero and future null. Every retained wall,
and every further acute wall, has nonpositive product with nu. The ideal
point therefore persists in every acute completion. Both retained pi/5
walls0,1 are incident to it.

A finite-volume completion would give a Euclidean cusp with a full
translation lattice. The product of those two incident reflections fixes
the orthogonal complement of their rotation plane and has linear trace
`d-2+2cos(2pi/5)`, which is irrational. A lattice-preserving linear map has
integer trace. This contradiction is the previously reviewed
[persistent-cusp obstruction](fresh_h9_fixed_core_cusp_obstruction_20261013.md),
also explained in the [weighted-binary review](fresh_h20_weighted_binary_20261014_dotted_completion_review.md).

The proof uses an arbitrary common strict future interior X, not a demand
that every extension retain a particular numerical center. It excludes
these recipes for all parameters and fields. It does not exclude the
remaining1119 nondegenerate conics or arbitrary right-S6 walls.

```sh
python3 scripts/fresh_nonarithmetic_h3_two_contact_20261016_affine_gate.py --verify
```
