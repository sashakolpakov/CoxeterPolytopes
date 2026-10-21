# Outside-author review of the whole-B bulk/endpoint theorem

20 October 2026. **PASS**. Full reviewed source:
`results/compact_whole_b_bulk_weak_packet_and_endpoint_pressure_20261020.md`,
SHA256 `e5a3c51858bb7b1e392b2828c874387be8d2044d5ae1be583a546b67eaddb54c`.
No correction is requested and no unresolved mathematical gap was found.

I did not coauthor this whole-B proof. I authored the imported A21
argument, so that dependency is disclosed rather than presented as an
unrelated outside invention. Here I independently rederived the
B-specific inventory, normalization, inequalities and final assembly.
No scientific job, worker import, optimization, Lean or TeX execution
was performed; accepted sources and the C28 packet were not changed.

## Complete bulk coverage

For 5<=i<=n-4 the retained A_(i-1) and B_(n-i) both have rank at
least four. A new node has at most one contact to each component.
The right B component already contains a 4 bond and has rank at
least four, so a finite extension can only attach 3 to its long end;
there is no omitted exceptional diagram of rank at least five.
The left possibilities are exactly the stated end 3/4, second-end 3
and proper E positions. A simultaneous completion of both components
has rank n and is the B path with simple endpoint joins.

I checked the zero-own-contact cases: the one full-path end 3 gives
a forbidden finite rank-(n+1) extension, internal 3 contacts have the
terminal-4 affine fork, end 4 has the two-terminal-4 affine path,
and two white contacts make the exact affine cycle. Zero contact
everywhere cannot be a unit spacelike outside normal. Thus a released
contact was not silently assumed positive.

The actual edge formula (1) follows directly from the retained white
inverse and C_ij=2 min(i,j). Own 4 with no right contact has the
terminal-4 affine path, and own 6 the exact affine 6,3 subpath.
For a near white-4/own-3 packet I independently obtained
`cosh L_i=(i-1)+i/sqrt2`; it is strong already at i=5, without an
A-type mate import. The other left/right cases meet the stated
38/25, 8/5 or 100/61 comparisons. Proper E positions start at i=6,
so their residual really is at most 1/4. Isolated caps use only the
local simple arms and D>=85/11 in the accepted A21 estimates; the
bulk B vertex supplies both arms before its terminal 4. This proves
that every weak bulk cap has precisely its own contact and helper
3 at node 1, including all native 5, >=7 and dotted possibilities.

## Native pair gap and actual mixed vertex

I rederived `A_i^2=2i gamma_i^2+2 gamma_i-1/2` and the uniform
gamma<6/5 seal. Weighted AM--GM with time weights b/a and a/b
gives exactly
`c_i.c_j>=1/2-ab+(a/b+b/a)/4>=1-ab>-11/25`.
The caps are distinct because their complete positive own contacts
are at different nodes. The actual wall-product gap therefore forces
orthogonality. The retained common-helper cost is
`t=(i-1)/(2i)`, and `I-tJ` has positive eigenvalues and an entrywise
positive inverse. Together with the actual old facets this gives an
actual finite vertex by positive acute projection.

The imported same-helper calculation uses only i>=5, A_i+A_j>4
and gamma>4/5; its decreasing squared upper bound gives
`sqrt(1+Dcal)/(W_i+W_j)<7/16` at all larger i. In B coordinates,
`sqrt(D_i D_j)-C_ij<1` and C_ij>=10 give the claimed extra loss
less than `W/(4sqrt10)`. Their sum is below .61. This excludes
adjacent weak indices without using an opposite-end assignment,
a partial Gram as a realization, or a finite-only label truncation.

## Weak-5 endpoint pressure

I reconstructed all four pressure vectors from the actual old-normal
relations. The finite-5 projection uses A1+A2 and its H4 completion;
7/8 use the negative-core point from the two orthogonal contacted
old walls; >=9 uses the single-contact positive projection; dotted
uses the actual dotted foot. Strict concavity and the positive
dotted coefficient gamma-gamma^-1 justify every substitution by
`gamma*g5+h1/2`, including the orientation of the resulting x1 bound.

The finite-label intervals give the asserted .8/.7/.65 bounds.
For the entire >=9 interval I checked A'>sqrt10, H_A>9/10,
H_delta<0 and the endpoint rational estimates, so this is not a
sampled endpoint argument. For dotted gamma in (1,6/5), the stated
derivatives make the difference strictly decreasing. Using precisely
the displayed rational radical bounds gives at 6/5 the strict lower
value

`1.41421*(4.03732-.552771-1.2*(3.162278-.61))`
`=.596535556834>.5965`.

The prefix minimizer at nodes 1,6,16 has all positive decreasing
slopes. Thus the lower constraints, rather than unjustified
equalities, give the energy bound. The strict seals for S6 and
4-sqrt6 check, and I independently expanded the displayed rational
sum to exactly `1.00043812051>1`. Replacing 16 by 17 increases its
last energy term: `(sqrt(t)-sqrt6)^2/(t-6)` increases with t.
Therefore index 5 is strong in every native case.

## Final whole-B exclusion and scope

At n>=26 both 21 and 22 are bulk. Nonadjacency supplies a strong
endpoint M=21 or 22, while 5 is strong. The reference prefix is
concave, with nonnegative dual coefficients, including its initial
linear part and constant continuation. I checked the strict initial
energy bound and the sum 35831/100000 of the sixteen reciprocal lower
bounds, giving E21>1.00971. At most eight internal omissions occur;
their segment interiors are disjoint. Since the square-root curvature
is less than 1/40, their total energy loss is strictly below
8/3200. Concavity survives and the remaining constraints still have
that interpolant as their minimum-energy vector. This validates the
final contradiction independently of the last four directions.

The result excludes a WHOLE B_n farthest vertex for n>=26. It is not
an arbitrary-spectator B-component theorem, an exclusion for n<=25,
or an unrestricted compact dimension ceiling. Finite spherical
classification and the pinned actual-point geometric lemmas remain
explicit inputs; no executable or formal geometry verification is
claimed by this written review.

Direct supporting sources read and rehashed:

| Source | SHA256 |
| --- | --- |
| `compact_whole_a21_farthest_vertex_exclusion_candidate_20261020.md` | `3abe11bb70ce677069279462a91111ed6b8ab8889e3ec9e656ac71c1094200fb` |
| `compact_a21_adjacent_weak_far_m_simple3_packet_20261020.md` | `1c4437c98d916d558081572dd2ecc12bc75eda1328a53633fe111120ec63dcb7` |
| `compact_farthest_vertex_mixed_replacement_inequalities_20261020.md` | `3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84` |

All supporting paths in this table are under `results/`.
