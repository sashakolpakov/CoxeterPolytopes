# Independent whole-A21 farthest-vertex proof review

20 October 2026. Reviewer: batch_guard. PAPER-ONLY mathematical review;
no scientific worker, optimizer, Lean, or TeX was executed.

**Verdict: PASS. No mathematical gap found in the stated theorem.**
A farthest vertex of a compact hyperbolic Coxeter21-polytope cannot have
WHOLE finite diagram A21. This is a geometric statement, without an
arithmeticity assumption. The review does not assert a dimension ceiling,
an extension to other ranks, or an extension to A21 with spectators.

## Bound source and dependencies

The complete source read and independently checked is
`compact_whole_a21_farthest_vertex_exclusion_candidate_20261020.md`,
SHA256 `3abe11bb70ce677069279462a91111ed6b8ab8889e3ec9e656ac71c1094200fb`.
Its earlier pending-review label is superseded by this review for these
exact bytes; the source itself was not edited.

The following imported sources were checked against the displayed pins:

* `compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`,
  SHA256 `3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`;
* `compact_an_projection_pressure_and_weak_cap_packing_review_20261020.md`,
  SHA256 `8b0a59493f1fa5dc7081e5e70beff2839ac5918bfcd50e089e05ae9b5546e43b`;
* `compact_a21_adjacent_weak_far_m_simple3_packet_20261020.md`,
  SHA256 `1c4437c98d916d558081572dd2ecc12bc75eda1328a53633fe111120ec63dcb7`;
* the preceding outside review
  `review_compact_a21_adjacent_weak_far_packet_batch_guard_20261020.md`,
  SHA256 `a64d90c7b01fda3de511c9230f8f85646e51701400e197d378e7daf341f58b42`.

All paths above are under `results/`. I contributed to the earlier general
mixed-projection mechanism and reviewed the adjacent packet. I did not
author this whole-A21 proof. Its new complete inventory, positive-dual
certificates, and final assembly were checked independently.

## Completeness and geometric interfaces

The unit-normal normalization, exact inverse A21 Gram entries, strict
concavity of the old dual heights, and energy
`(1/2) sum (h_(i+1)-h_i)^2 < 1` are consistent. Both positive-definite
face projections and negative-delta positive combinations give actual
feasible hyperbolic points. The dotted-pair foot has the asserted norm
and all-wall signs. There is no replacement of full tangents by half
tangents in the pressures used here.

For every central direction 5 through 17, the retained paths both have
rank at least four. The finite endpoint classification covers isolated
walls; one white3, including the D/E internal possibilities; one white4;
and the two-white3 bridge. The rank21 bridge must be A21. The separate
zero-released-contact discussion excludes every branch, rather than
assuming the released product is positive. The affine/nonpositive
subdiagrams used later are retained at the relevant actual endpoint.

The isolated finite and dotted cases all give the stated strong floor.
In particular, the finite-label-at-least8 estimate uses the FULL function
H: if b=sqrt(D-1) and t=tan(theta), then
`H/gamma = b sqrt((b-t)/(b+t)) > b-t`.
It is not obtained from the weaker negative-delta bound. The isolated7,
isolated5 H4, and dotted loss estimates also check. The near4 cases have
actual strong mates; their separation from an opposite strong direction
gives the strict contradiction claimed in Section 3. The remaining weak
directions are precisely covered by the imported adjacent packet, so
they cannot be adjacent.

## Positive duals and final energy

I checked every helper row for labels5,7,8, finite labels at least9, and
dotted products, including the separate left/right helper locations and
the strong16/strong17 alternatives. The projected variances, height
demands, monotonicity estimates, and signs of all old constraint
multipliers agree. Thus these are positive-dual lower bounds, not
equality-constrained minimizations with unverified multiplier signs.
The smallest displayed certificate gaps remain positive: the finite8
gap is 9/40000 and the finite-at-least9 gap is 219/39062500.
Consequently indices5 and17 are strong.

The radical boxes used for the final scalar seal can each be checked
by squaring positive rational endpoints. They give the direct lower sum
`0.941519618 + 0.06567285 = 1.007192468` for the full central-floor
energy. For q(t)=sqrt(t(22-t)/11), the exact identity is
`q''(t)=-q(t)^(-3)`; hence the asserted curvature bound and isolated-hole
loss of less than 1/1600 follow. Remaining weak indices are nonadjacent
and number at most six. Their affected interpolation intervals have
disjoint interiors, their losses add, and the surviving concave envelope
retains nonnegative old multipliers. The resulting energy is strictly
greater than `4013/4000 > 1`, contradicting the farthest vector norm.

This completes the outside mathematical audit of the whole theorem.
No changes to its source or imported evidence were needed.
