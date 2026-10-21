# Independent scalar and composition review of dimension69

The [dimension69 proof](fresh_ideal5_weighted_dimension69_20261017.md) passes
the independent exact reconstruction and mathematical review described here.
No remaining gap was identified in the reviewed scope. The
[review source](../scripts/fresh_ideal5_weighted_dimension69_20261017_review.py)
uses only the Python standard library and imports no project arithmetic.
Its [saved review certificate](fresh_ideal5_weighted_dimension69_20261017_review.json)
binds the frozen primary and its direct inputs. Final generation and a
separate fresh saved replay both passed; their scope is recorded below.

## Independent scalar calculation

At ambient dimension70 the cusp rank is69. The reviewer enumerates profiles
by choosing an optional high rank, an optional low rank and an optional12
first, then filling the remaining rank with copies of11,13,14,15,16,17.
The primary uses its inherited repeated-ranks-first enumeration. Exact set
equality gives the same1,339 terminal profiles, without duplicates.

For each profile the reviewer computes the polynomial
`H(t)=product_i(1+t+...+t^r_i)` by sliding-window convolution. It checks
symmetry, the vertex count V, and independently recomputes G by the
cumulative bounded-composition identity

    G=floor(d/2)*V-2*sum_shift a_shift*binom(floor(d/2)-shift+s,s+1),

where `sum a_shift*t^shift=product_i(1-t^(r_i+1))` and s is the number of
factors. It also checks the absolute-deviation identity

    curvature/V=3d/4-(5/2)*E|sum U_i-d/2|.

All coefficients, vertex counts, G values, normalized curvature, costs and
margins agree exactly with every primary row. The minimum margin is
**1819/51840**, attained at `(11,11,11,11,11,14)`. The finite slack is
**0**, and the strict global additive constant is345. The ambient69
boundary profile `(11,11,11,11,11,13)` has margin **−2029/13608**; it is
a formal obstruction to this comparison, not a claimed realized cusp.

Only the rank13 and14 coefficients change, becoming43/7 and32/5. The
complete concentration domain remains valid: low-rank concentration has
unchanged linear cost, two12 ranks become11 and13 with exact cost gain
19/182, and the constant-cost high ranks concentrate to at most one rank
above17. Every nontrivial move increases the bounded integer squared-rank
sum. No further concentration or elimination of ranks13 through17 is
assumed. The reviewer also reconstructs the unchanged fifteen all-rank
affine cost branches and their endpoint/limit maxima, rather than inferring
an all-rank maximum from a finite sample.

## Independent affected-map and projection checks

The [new local theorem](fresh_b14_b15_singleton_completion_transfer_20261017.md)
deletes one pair from each of the order14 and15 affine B rules. The reviewer
independently builds both canonical target graphs, their distances, feature
sets and two automorphisms. It verifies that each changed pair is an
invariant distance10 singleton, and reconstructs the paid sets of sizes86
and96. All elementary triangle, short8 and basic ideal pairs remain units.

To determine the complete affected support, the reviewer derives marked
images directly from all506 saved order14 maps and548 saved order15 maps
in the [frozen complete small-factor certificate](fresh_distance8_ideal5_small_factor_20261017.md).
These1,054 maps occupy422 target/context rows. The derived marked-image
sets agree with the saved sets, and select precisely46 and18 affected maps
in23 and11 contexts. Five contexts occur for both targets: there are29
distinct literal contexts and34 target/context occurrences. Every affected
map is checked injectively against the labelled old core Gram and appears
exactly once in the new local certificate.

For the thirteen distinct old core Grams, independent Fraction Gaussian
elimination reconstructs the full inverses. Exact congruence elimination,
including a two-by-two hyperbolic pivot when necessary, checks their Lorentz
signatures. For each target/context occurrence the reviewer checks the
supporting-white positive principal and all four strict edge principals,
reconstructs the preferred finite corner and its connected components,
and derives the complete necessary contact domain from the finite-block
inverse diagonals and the two omitted blacks. The result agrees with the
local certificate and the original overlap record.

The reviewer then checks that the two end sets partition this domain,
their full inverse blocks are entrywise positive, and touching both ends
supplies a retained short path. It derives the shortcut length from the
actual old core, rather than accepting the saved bound. The other old
corner components are orthogonal to the marked component, and the two
omitted blacks are absent from the preferred corner.

Every affected map forces the target's missing doubled terminal to touch
an old marked node i. Because the terminal is outside the displayed
component but adjacent to it, it cannot be another old corner component;
an omitted black also cannot belong to the actual corner. Thus it is a
new supporting white. In the norm-two normalization its required contact
has square2. All64 exact calculations give `2(H^-1)_ii>2`, with minimum
5/2, and agree with the local certificate's witnesses.

The mathematical argument permits arbitrary further contacts. Contacts
on both ends give the already-covered short path. If all contacts lie in
one end, its entrywise-positive inverse block makes the forced term a
lower bound for the full projection square. The core is Lorentzian, so
the perpendicular residual has nonnegative square and a norm-two root's
projection cannot exceed2. This contradiction excludes the affine
occurrence. No unspecified Gram entries are assigned zero, no indirect
tail cutoff is imposed, and no displayed finite subdiagram is silently
treated as a whole component.

The complete old core and induced-map enumerations remain frozen inputs.
This review checks all saved target-map images and every affected literal
projection; it does not independently repeat those original exhaustive
censuses.

## Global composition and zero finite slack

The full local domain split is unchanged. Triangles and short8 paths retain
their required units. Low-core alternatives have whole component order
at most12; the high-disjoint branch is excluded; protected/basic ideal
payments have distance at most9; finite and universal cycle alternatives
retain their rules; and the universal immediate whole-affine alternatives
have order at most13. The complete remaining overlap support is the one
checked above. Thus triangle and quadrilateral totals remain at least2
and1, and no finite coupon or extra finite cost is introduced.

The [finite budget](fresh_finite8_d12_b10_coupon_budget_20261017.md) remains8n
through dimension75. Its previously independently checked13,034 endpoints
are inherited through the [dimension70 review](fresh_ideal5_weighted_dimension70_20261017_review.md)
without repeating that unchanged enumeration. The new order14 maximum
is86; the order15 maximum is96, shared by B and D. Every other small-order
cost and the all-rank tail are unchanged. The actual cusp factor-fiber
count multiplies a factor's pair weight by V/(r+1), including at nonsimple
cusps.

At n=70 the shape inequality has finite coefficient
`beta_70=70*64/8=560=8n`. Subtracting the face-payment upper bound from
the [shape lower bound](fresh_cusp_link_shape_curvature_20261017.md) leaves
zero times the finite-vertex count, strictly positive coefficients at all
ideal vertices, and the unconditional additive term345. The difference
is therefore positive even with no ideal vertices. Positive finite slack
is not required. This excludes dimension70, and the completed dimension70
theorem excludes every higher dimension, establishing the unrestricted
ceiling69.

No arithmeticity, simple-cusp, factor-count or facet-count hypothesis is
added. The construction problem and compact-arithmetic models are separate;
this review supplies no new example, external refereeing or formalization.

## Verification status

Independent scalar preparation, complete primary/local comparison, final
hash-bound generation and separate fresh saved replay all terminated
exit0/PASS. Syntax, all27 direct source bindings and all8 local report
links pass. The full primary source and written proof,
the new local source and written argument, and the global composition
received a complete mathematical read; no correction was required. The
parent independently read this reviewer's full source and mathematical
reconstruction, including scalar enumeration and all three curvature
identities, inverse/inertia and actual-corner routines, contact-domain and
map coverage, and the zero-finite-slack composition. No gap was identified;
that read did not repeat the reviewer's1,339-profile calculation.

The bound primary hashes are source
`013692566433013fa44efad68c9abb21c7bdfa666f56a3e51e8da4697cdacf0d`,
data `70a9bc7d0d5e13252b6e9e7743fb989759bc5d4521901d548b5130bda39ecb8c`,
and report `600a9acb16136c1e30a1d591a1ef6375b7c841d1378968d3aa25b7bf72be5c4e`.
All validation is internal. The primary and original censuses are frozen
dependencies with the distinct verification scopes described above.
