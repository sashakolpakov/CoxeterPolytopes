# Outside audit of the whole-A21--A29 plateau transfer

20 October 2026. Reviewer: resource_cleanup. **PASS.** This is a full
written mathematical review, with no scientific execution, imports,
optimizer, build, or formal verification. The source was not changed.

## Source binding and scope

Reviewed in full:
`results/compact_whole_an_geometric_a21_plateau_transfer_20261020.md`,
SHA256 `8a58b7747e1a3d8c9affe7df17038781a62983ddd70ccb670502c6f40b1f9cee`.

I also read its unchanged whole-A21 proof
`compact_whole_a21_farthest_vertex_exclusion_candidate_20261020.md`,
SHA256 `3abe11bb70ce677069279462a91111ed6b8ab8889e3ec9e656ac71c1094200fb`,
the outside A21 review `618e255766d77a00ea4457bf3f5b6d983ebb845bcd8352a259d291503dbf7c3f`,
the general packing/mate source `8b0a59493f1fa5dc7081e5e70beff2839ac5918bfcd50e089e05ae9b5546e43b`,
and the old adjacent packet `1c4437c98d916d558081572dd2ecc12bc75eda1328a53633fe111120ec63dcb7`.
These paths and hashes agree with the stated inputs. I co-developed the
earlier projection-pressure/packing interfaces; I did not author this
plateau transfer. Its new bounds and exact transfer were checked anew.

The result concerns a WHOLE A_n diagram at a farthest vertex, for
21<=n<=29. The extension to the phrase all n>=21 uses the explicitly
imported universal compact ceiling29. No spectators, arithmeticity, or
new ambient dimension ceiling are asserted.

## Uniform inventory and simultaneous packet check

The two retained paths still have rank at least4, and the two-path
completion has rank at least21. Thus the accepted finite-white and
zero-own-contact proofs cover every central endpoint. I checked that
the isolated estimates need only the lower dual norm. The new scalar
comparison clears to the printed concave quadratic; its endpoint
values on [5/2,31/8] are positive.

The near4 ratio increases with its side parameter m. At m8 its value
is `2sqrt(2)/3+7/12>761/500`. At m7 and N22/23 the dual norm falls
in the old [3,25/8] comparison interval. This leaves exactly the stated
mate locations, including index6 when N>=24. The opposite-anchor
correlation is <1/3 for N22/23 and <=6/(N-6)<=1/3 thereafter, which
suffices for the strict normalized-height contradiction.

For literal packets I rederived the weighted AM--GM formulas(4), with
both helper orientations. They give product0 for a common helper and
only0 or -1/2 for opposite helpers. In particular the new lower bound
>-7/10 already excludes the native4 product. The dotted upper cutoff
is valid: after clearing the positive denominator its difference from
`4g-2` is `[2m(g-1)-1]^2`.

The outward zero-product determinant clears to
`(2/5)C-24/25+8/(5N)>0`. I independently expanded the reversed
determinant and obtained precisely equation(5), including its constant
`(N+2)/4`. Its derivatives are positive on the stated rectangle;
at the boundary pair(4/5,9/10), concavity in the integer position gives
the lower value `(509N-9974)/500>=306/125`. The remaining two5 case
is exactly covered by the geometric packing input. Thus the mixed
Schur matrices are positive definite with positive inverses in every
remaining branch.

The mixed-point bounds `7/16` and `3/10` remain valid with only the
displayed side-length and time lower bounds. Rationalizing the new
E0-minus-norm loss gives <3W/32. The resulting weighted losses
<17/32 or <63/160 exclude adjacent weak literal packets. The ensuing
near4 removal is therefore not circular.

## Exact positive-dual transfer and final energy

In potential coordinates, inserting a plateau changes no nonzero
first difference. It is a linear isometric lift. More explicitly,
for i<s the potential of `omega_i-(i/s)omega_s` equals
`2min(i,j)-2ij/s` on j<=s and0 on j>=s, independently of N.
This verifies the claimed end-supported pressure residuals, their
unchanged variances and projections, for both helper assignments and
both strong-right-anchor alternatives. It avoids falsely identifying
a lifted individual omega with the larger-rank omega.

The old peak coefficient splits into two nonnegative slope drops at
11 and N-11; their sum is the original coefficient. All corrected
end-anchor coefficients retain the accepted positive signs. The
actual time formulas have the stated negative constant-over-N terms,
so the selected pressure targets increase with N, while their local
Schur deltas stay fixed. Consequently every old positive-dual weak5
certificate transfers with its exact norm and at least its old margin.
Concavity and nonadjacency supply both peak-anchor floors, even when
one of those indices is weak. Reversal excludes weak(N-5).

For N>=23 the two actual peak heights exceed or equal the reference
`sqrt(11)-.61`: a weak peak has strong neighbors whose dual squares
are at least `260/23>11`. The lifted full reference has unchanged
energy >1007/1000. Only the five sites6,...,10 in each flank can be
holes, and nonadjacency permits at most three per flank. Their old
interpolation intervals have disjoint interiors, including at plateau
endpoints; each has the unchanged loss <1/1600. The surviving curve
is concave with nonnegative dual coefficients. Thus its lower energy
is >4013/4000>1, the required contradiction.

No mathematical gap or source correction was found in the claimed
rank21--29 transfer. No executable or formal-geometry coverage is
claimed by this review.
