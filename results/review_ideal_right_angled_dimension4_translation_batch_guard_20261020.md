# Independent full manuscript and model-source review

20 October2026. **Written mathematical translation review PASS.**
No scientific execution, model replay, Lean compilation or TeX build
was run by this reviewer. The following three sources were read in full:

|Source|SHA256|
|---|---|
|`manuscripts/ideal_right_angled_dimension4.tex`|`c780ba2f7243553b0deaa89bece75964f9bb2903888536c88b501e698aa06db9`|
|`results/ideal_right_angled_dimension4_proof_20261020.md`|`9e6f7eac5438c3a71c8c90bc5c615dcdac51072dfb7869d9654ef3c76e9d8768`|
|`scripts/check_ideal_right_angled_dimension4_model_20261020.py`|`bb506ec1ee98dee1429d6e668553af40316ebdeda6b59b82bbf03b8f3a6dc261`|

This review supplements the unchanged complete adversarial review
`review_ideal_right_angled_maximum_height_kissing_batch_guard_20261020.md`,
SHA256 `053fff76e71ea51a0be89bc2938b59e3504142a8ae697d0a303e37bf0f9e5550`,
which binds the original95ef candidate and independently checks its
published inputs. None of those earlier bytes was edited.

## Complete mathematical translation

The authoritative note and manuscript retain the original complete
argument. The journal reference now correctly names Kolpakov's Theorem4,
printed p1945, with its metric equality case. Musin's input is precisely
the finite60-degree-code cardinality in S3, without a code-uniqueness
assumption. Vinberg's all-wall sign is explicitly invoked for nonadjacent
facets, rather than inferred only from the actual ridge angles locally.

I rechecked the following directly in the TeX, rather than relying on
the prose proof note: unit-normal height normalization; existence of a
positive finite maximum; q=v/H0+n with q^2=-1 and q.v=-H0; feasibility
against every wall; the exact intrinsic height ratio at actual orthogonal
neighbors; injectivity of their radial directions; both inequalities
in the S3 code bound; ideal finite-volume right-angled inheritance by
the four-facet; the24-facet cardinality squeeze; and the explicit opposite
normal-pair proof that q is the center. All are correct.

Height-zero neighbors remain allowed until the center proof excludes
them. The equality argument does not incorrectly apply a uniqueness
theorem to the kissing code. The repeated local argument applies to
every globally maximal-height facet, so its set is closed under actual
ridge adjacency. The explanation of connectedness via the compact
projective closure correctly distinguishes those graph edges from
cusp-only contacts. Propagation reaches an actual height-zero cusp
facet and contradicts H0>0.

The higher-dimensional tail uses actual five-face inheritance of entire
ideality, finite volume and right angles. It does not infer a universal
decomposition or impose an arithmetic or facet-spectrum restriction.
The24-cell supplies the sharp four-dimensional example. No theorem
in the manuscript is stronger than the reviewed argument supports.

## The two figure panels

The first panel has upper half-plane circles of radii1 and2, centered
at0 and sqrt5 on the boundary. Their center-distance square is5=1+4,
so their hyperbolic walls meet orthogonally. The displayed q=(0,1) is
the top of the smaller circle and lies on the exterior side of the
larger. The distance to the larger wall has

    sinh distance=(5+1-4)/(2*2*1)=1/2.

Thus the caption and its height ratio H'=H0/2 are exact. The dashed
upward line represents the chosen ideal endpoint, not another facet.
The second panel is explicitly combinatorial propagation along a facet
path. Its placement of graph nodes is not used as a hyperbolic metric
claim. Both captions accurately state their scope.

## Exact model checker: source-only scope

The checker constructs24 distinct unit D4-wall normals and checks all
276 unordered Gram entries, twelve opposite pairs,24 displayed feasible
null vertices with six incident walls and three opposite cusp pairs,
and an actual strict facet point for every wall. Its finite facet point
(2,r) has square-2, wall product0, and strictly negative products with
all other walls. These checks agree with the written standard model.

Its scalar identity

    (1+A)(1+B)-4AB=(1-A)(1+B)+2A(1-B)

has correct coefficients and supports the second code inequality for
A,B in[0,1]. The script does not claim that an execution of this identity
proves all real domain facts. It explicitly disclaims universal geometric
verification and completeness of the displayed vertex enumeration.
That is the appropriate boundary: the universal proof is the written
argument and its named published inputs. I found no source flaw, but
have not executed the checker or asserted its runtime result.

## Boundary of this review

This is a full source-bound mathematical review, including both figure
captions. It is not a TeX layout/visual inspection, an execution ledger,
an independent proof of Musin's theorem, a reproof of Kolpakov's theorem,
or a novelty/priority search. Conditional Lean logic has a separate
source and must receive its own resource-gated compilation and audit.
