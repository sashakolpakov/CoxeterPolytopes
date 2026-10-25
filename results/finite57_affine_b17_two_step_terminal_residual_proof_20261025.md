# A missing doubled terminal after two white roots

25 October2026. Argument by `/root`, prepared for independent review.
This is a local proof proposed for the existing finite57 verification.
No new dimension bound, computation, application or example is claimed.

## The precise setting

Use norm-two outward normals in a real Lorentz space with one negative
direction. The roots of an original quadrilateral are divided into
supporting whites and four boundary blacks. Write the opposite black
pairs as `(a,c)` and `(b,d)`. Thus `a,d` and `c,d` are adjacent pairs.
All actual facet products are nonpositive. A nonzero product between
distinct norm-two Coxeter roots has magnitude at least one.

Suppose a displayed old core has seventeen roots, nonsingular Gram matrix
`G`, and exactly one negative direction. Let `J=G^{-1}`. Among these old
roots, those incident to the selected corner are precisely the fifteen
other than `c,d`. Assume an actual whole affine B factor of order17 completes their
given full native map. Use the canonical chain `0,...,14`, with a doubled
bond `13--14` and two simple leaves `15,16` at node0. The given map omits
exactly `13,14`, sends `a` to9 and `b` to one of15,16, and retains the
complete native edges and nonedges. Let `r` be the old root at node12.
The old roots at10 and11 are also present.

Let `u,t` be the roots at the missing nodes13,14. They are genuine new
supporting whites. Indeed, they are different roots of the actual whole
factor, while the other original blacks `c,d` do not contain this selected
corner. In particular they cannot stand in for the missing white roots.
The native map gives

* `(u,r)=-1`; `u` is orthogonal to every other old selected-corner root;
* `t` is orthogonal to every old selected-corner root;
* `(u,t)=-sqrt(2)`.

Assume the following exact inverse inequalities, in this physical unit
gauge:

    J_cc > 2,    J_dd > 2,    J_rr > 1,    J_rc >= 0.

These are local numerical hypotheses. They must be checked or bound to
the exact original inverse; a core name or matching count is insufficient.
There is no assumption that `u` has no contact with `c`.

We prove that this completion is impossible if the same actual face has
no adjacent-black path of length at most8 with white interior. The
independent current S8 theorem supplies this hypothesis on an unpaid face.

## Projection and the terminal root

The old span is nondegenerate and has one negative direction. Its
orthogonal complement in the ambient Lorentz space is therefore positive
definite. For a new norm-two root `v`, let `q` be the column of its old
products. Its projection onto the old span has square `q^T J q`, and its
orthogonal remainder has square `2-q^T J q >= 0`.

The only possible old contacts of `t` are with `c,d`, since the full native
map prescribes zero products with all fifteen other old roots. Contact
with both would give the actual path `c--t--d`, with white interior and
length2. This is excluded by the no-short-path hypothesis.
If just `c` is contacted, its product magnitude is at least one, so the
projection square is at least `J_cc>2`, impossible. The same argument uses
`J_dd>2` when only `d` is contacted. Consequently `t` is orthogonal to the
entire old span, and its square in the positive definite complement is2.

## The preceding root and the contradiction

The root `u` cannot contact `d`. Such a contact would produce the path

    a(9) -- old(10) -- old(11) -- r(12) -- u(13) -- d,

of length5. Its interior consists entirely of actual whites: the only
selected blacks have images9 and15/16, and `c,d` are outside the fifteen
mapped old roots. Its endpoints `a,d` are adjacent original blacks. This
is precisely a forbidden same-face short path, not an abstract graph
shortcut involving an internal black.

Write `-(u,c)=s>=0`. The full old product column of `u` thus has entries
`-1` at `r`, `-s` at `c`, and zero elsewhere. Its projection square is

    J_rr + 2 s J_rc + s^2 J_cc >= J_rr > 1.

If this already exceeds2 there is an immediate contradiction. Otherwise
the remainder of `u` in the positive definite complement has square
strictly less than1. As `t` is entirely in that complement, Cauchy's
inequality gives

    (u,t)^2 <= 2 (2 - projection_square(u)) < 2.

The required doubled bond has `(u,t)^2=2`, a contradiction. This proves
the stated exclusion while allowing every possible old contact of `u`
and every ambient complement dimension.

## Application boundary

The failed selected515 map is a candidate for this lemma. Reading its
stored coefficients gives `J_cc=7/2`, `J_dd=4`, `J_rr=3/2`, and
`J_rc=3/2` for the corresponding physical nodes `c=1,d=12,r=4`.
This was retrieval of existing data, not a new inverse computation.
Every actual application must still bind its original opposite pairs,
physical gauge, full native map, old nonsingular Lorentz span and exact
inverse. The analogous other marked leaf is covered only after checking
the same hypotheses; its applicability is not inferred from a high index.

A pure application adapter may verify these hypotheses for a given map.
It does not assert that a short path was observed. It proves an exclusion
under no short path; composing with the independently accepted S8 theorem
then rules out that completion in an unpaid face. Source review and root
acceptance must precede any such execution or payment credit.
