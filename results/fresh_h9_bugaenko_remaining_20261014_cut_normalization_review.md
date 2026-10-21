# Review of nested two-cut normalization

**PASS.** Source/proof review and a fresh saved replay of
[the checker](../scripts/fresh_h9_bugaenko_remaining_20261014_cut_normalization.py)
passed all 6,160 literal vertex cases. The [source/result bindings](fresh_h9_bugaenko_remaining_20261014_cut_normalization_review.json)
record the reviewed versions. No independent full polytope enumeration
is claimed or required for this local identity.

A truncation is stellar subdivision of its facet support in the dual
nerve. If the first cut has support `S`, it makes `S` a minimal nonface
while preserving every proper subset. A second subdivision with support
`T` can destroy a proper subset of `S` only if `T` lies in that subset.
Thus, unless `T` is a proper subset of `S`, the large minimal nonface `S`
survives. This includes supports containing the newly introduced facet.

For `T` properly contained in `S`, the checked identity exchanges the
order: first truncate `T`, introducing `b`, then truncate
`U={b} union (S minus T)`, introducing `a`. The second support is a face
after the first normalized cut: every old vertex containing `S` produces
a new vertex containing `U`. The displayed vertex identity handles all
three cases `T` absent, `T` present but `S` absent, and `S` present.

The normalized first missing face `T` survives because the second support
contains `b` and is not a subset of `T`. The normalized second cut creates
the minimal nonface `U`, retaining every proper subset. Consequently any
compact hyperbolic Coxeter realization forces both `|T|<=5` and `|U|<=5`
by the complete Lanner bound. If the second original support itself had
size above 5, its newly created minimal nonface would already exclude
compactness. If the first original support has size at most 5, the
original order is already among the small-cut cases.

These cases justify normalization of every potentially admissible pair
of proper face truncations of codimensions 2 through 9 to two cuts of
codimensions 2 through 5. The extension of the particular Bugaenko
wedge/two-cut exclusion still depends on the separately certified
exhaustive small-cut catalogue and its obstructions. This local review
does not by itself certify those additional searches.
