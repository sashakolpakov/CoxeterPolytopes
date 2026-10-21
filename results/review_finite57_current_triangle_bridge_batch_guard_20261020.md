# Current-rule triangle bridge: source audit and proposed insertion

This is an internal mathematical source audit, not external refereeing or a fresh computational replay. I reconstructed the contact obstructions and compared the current marked rules with the selected triangle witnesses. Earlier work on this manuscript and exchanges with the parent agent are disclosed; no prior PASS was used as a mathematical premise.

## Bound sources and result

- Canonical source inspected, unchanged: `manuscripts/finite_volume57.tex`, SHA-256 `89bd7658a909c566185544484c4a95f8d176fb5fe476cffa7befad81602ca043`.
- Proposed supporting fragment: `results/finite57_current_triangle_bridge_insert_20261020.tex`, SHA-256 `09c13141ada8bd23373af6da39e18a3a8b3b626701d57b6dce57071ca3960348`.

**Finding:** Appendix `app:routing`, “Triangles and their four exceptional markings” (currently lines 1031–1074), needs an explicit current-rule bridge. The inherited claim that a path of length at most five pays a unit has exceptions in affine cycles of orders ten and eleven. The order-ten affine-D opposite-fork reduction also has distance seven. These pairs cannot occur at a corner of an actual triangle. Thus this is a missing written justification, not a counterexample to triangle payment two. The supporting fragment supplies that justification and the preservation comparison. Insert it immediately after the subsection heading, before the historical minimal-core census. It is not a new paper and has not been inserted into the canonical source.

## Independent mathematical checks

Normalize all actual normals to square two. At either other corner of the triangle, the deleted affine factor together with the third side is a principal positive-semidefinite submatrix of the actual corner diagram. If its finite deleted-factor Gram matrix is (B), its nonnegative contact vector (p) satisfies (p^TB^{-1}p\leq2). Entrywise nonnegativity of (B^{-1}), and the native lower bound one on a nonzero contact, exclude every contact position with inverse diagonal greater than two. This principal-submatrix argument allows arbitrary other supporting roots and cusp factors.

For cycle orders ten and eleven, with marks (0,5), the two allowed sets are respectively

- order ten: `{1,2,8,9}` and `{3,4,6,7}`;
- order eleven: `{1,2,9,10}` and `{3,4,6,7}`.

These follow from (A_{m-1}^{-1}(k,k)=k(m-k)/m), whose first forbidden diagonal is (21/10) or (24/11). The other deletion also forbids the removed mark, by diagonal (5/2) or (30/11).

For order-ten affine D, use spine `0..5`, leaves `6,7` at `0`, and `8,9` at `5`. For marks `6,8`, deletion gives finite D9. The only allowed sets are `{0,7}` and `{5,9}`; the removed marks have inverse diagonal (9/4). All three cases force the third side orthogonal to the entire varying affine factor. Its positive affine-kernel combination is a nonzero future null vector representing the actual ideal corner. The third side is nonincident to that corner and therefore has strictly positive height there, contradicting orthogonality. This uses the actual triangle, not an abstract labelled diagram. The relevant symmetries cover every affected pair.

The other printed current affine B/C/D replacements of orders ten through seventeen have distance at least eight. The larger-order selected length-six/seven paths retain their branch or nonsimple feature within five steps of a mark; they cannot lie in an affine-A cycle. Exceptional affine factors have order at most nine and retain all units.

The four alternate pairs in the displayed ten-node tree remain protected. For three markings, the opposite affine-E8 kernel forbids additional supporting contacts to its nodes. For the fourth, the finite-E8 inverse has diagonal `(30,8,14,4,20,12,6,2)`, so additional contacts can occur only at the last node. Thus the fixed simple leaf and its degree-two neighbor persist. No classical affine diagram has this feature. In a finite component the protected leaf lies at distance `2,1,0,0` from a mark in the four cases, preventing a deep finite deletion. The other selected finite six/seven paths likewise retain a feature within five steps; the sparse rank-12..17 finite replacements have longer distance. No selected triangle unit is lost under the current changes examined.

## Actual proof homes consulted

- `results/fresh_order10_cycle_half_transfer_20261018.md`, “Triangles and simultaneous priorities,” and its same-stem Python source: order-ten allowed-contact obstruction.
- `results/fresh_cycle11_half_pair_transfer_20261018.md`, “No triangle has an affected corner”: order-eleven obstruction and positive-height argument.
- `results/fresh_order10_d_opposite_fork_fractional_transfer_20261018.md`, “Triangles, current priorities and scope,” and its same-stem Python source: D9 inverse diagonals and opposite-fork sets.
- `results/fresh_finite7_ideal5_triangle_payments_20261017.md` and `results/fresh_ideal5_finite8_triangle_payments_20261017.md`, with their same-stem Python sources: the historical 1025 marked-core census, selected features, and four E8 alternate witnesses.
- `results/fresh_order10_b_terminal_half_transfer_20261018.md`, `results/fresh_order11_bcd_half_pair_transfer_20261018.md`, and `results/fresh_order11_bc_zero_transfer_20261018.md`: remaining modified distances and short-unit preservation; compared with current Appendix `app:weights`.

The historical finite core census and its reported outputs were read, not rerun. This audit supplies the mathematical interface from those selected outputs to the current weights; it does not independently regenerate the census or audit all quadrilateral, forest-maximization, or higher-tail computations. No worker, project import, Lean, TeX, solver, or scientific job was run. Canonical sources and frozen evidence were not edited.
