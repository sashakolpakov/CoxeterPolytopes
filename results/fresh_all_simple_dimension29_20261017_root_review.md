# Root review of the mixed-simple dimension-29 composition

**PASS — internal mathematical and source review.** The
[completed proof](fresh_all_simple_dimension29_20261017.md) excludes dimensions
n>=30 for finite-sided, finite-volume hyperbolic Coxeter polytopes simple
at every vertex. Simplicity includes every ideal vertex. No arithmetic
assumption is used. This review does not apply the conclusion to nonsimple
cusps, claim novelty, or constitute external refereeing or formalization.

The reviewed composer source has SHA256
`6a7c81654008c5fdcd92948b9ab4eb9ad72abf25dff9db3735ad0b9f821ce691`;
its data has SHA256
`28ab979d7daf3c59e9f8f8b610d402e8a0f655ba0707b3c0a992b8db5a9d6927`;
its final proof report has SHA256
`f20903a8a54f98d8ac4600747c091e4f85f113acd108dfaca55b014d766fdf61`.

## Geometry and the case partition

For an all-simple finite-volume polytope, each ideal vertex has one connected
affine diagram on n nodes. At n>=30, exceptional affine diagrams are too
small, leaving only classical cycles and trees. A quadrilateral therefore
either has four finite corners, has an affine-cycle corner, or has an
affine-tree corner. These are exhaustive branches, and every cited lemma
requires dimension at most18. The finite and ideal rules use the same marked
corner convention and nonnegative weights. The vertex budgets are at most3n;
the face payments are at least two for triangles and one for quadrilaterals.
For k>=5, nonnegativity supplies the required lower bound5-k.

I read the corrected [affine-tree proof](fresh_simple_cusp_affine_tree_completion_20261017.md)
and source, including the all-rank Green identity, the four terminal
replacement profiles, the D2 outer component, and the D6/D7 finite extensions.
The corrected fork-leaf exceptional opposite components are finite D8 and
finite E8, with the bond4 alternative containing a proper affine obstruction.
The separate exceptional Schur calculations use those correct descriptions.
The [terminal/E8 audit](fresh_simple_cusp_affine_tree_terminal_audit_20261017.md)
has its own physical-matrix arithmetic scope. This root review does not claim
to duplicate all82,542 producer formula probes.

I also compared the [finite-corner transfer](fresh_simple_finite_quad_vinberg_transfer_20261017.md)
with Vinberg's original §6.5–11. Its geometric replacements pass: the actual
compact face gives indefinite opposite-pair blocks and a full star of corank
one; orthogonality to an actual indefinite block gives positive definiteness;
the private blocks and their complements have the stated size bounds.
Lemma7 supplies connectedness before bounding an elliptic high-bond block
by four nodes. The enlarged neighbor blocks have at most ten nodes; larger
orthogonal end pieces use the geometric replacement without a size cutoff.
The dependence coefficients on the private high-bond block prove that the
quadratic determinant coefficient is nonzero. The small cyclic-core affine
tests in §6.11 are distinguished from the later complete finite-corner
classification. No order-n affine diagram is discarded merely for being
affine. Published local tables remain explicit inputs, as in the separate
[source review](fresh_simple_finite_quad_vinberg_transfer_20261017_review.md).

## Independent final inequality

Let f0 and f2 count vertices and two-faces, and let h_i be the h-vector of
the dual simplicial polytope. The projective closure is a simple convex
polytope, so the standard identities and inequalities apply:

`f0=sum h_i`, `f2=sum binom(i,2) h_i`,
`h_i=h_(n-i)>=0`, and `h_0=h_n=1`.

Every vertex belongs to binom(n,2) two-faces. Thus the sum of the required
face payments is `D5=5 f2-binom(n,2) f0`. Pairing the coefficient for i
with that for n-i gives the following independent derivation:

* For even n, its minimum is `beta=n(n-6)/8`, and the coefficient excess
  is `(5/2)(i-n/2)^2`.
* For odd n, its minimum is `beta=(n-1)(n-5)/8`, and the coefficient excess
  is `(5/2)((i-n/2)^2-1/4)`.

These excesses are nonnegative. At both endpoints they are positive, and
the corresponding h entries are one. Consequently `D5>beta f0` is strict.
In particular strictness is retained at n=30, where beta=3n.

Substituting `n=30+2t` gives `8(beta-3n)=4t^2+60t`;
substituting `n=31+2t` gives `8(beta-3n)=4t^2+64t+36`.
These are nonnegative for every integer t>=0, without a dimension cutoff.
The common weight rule would therefore give

`3n f0 >= total payment >= D5 > beta f0 >= 3n f0`.

This contradiction proves the stated dimension bound29. At n=29 the scalar
inequality alone does not contradict the budget, and the proof claims no
lower threshold. No unresolved step was identified in this composition.

## Separate construction review

A narrow source and conceptual read of the
[central-truncation obstruction](fresh_nonarithmetic_central_truncation_domain_20261017.md)
also passes. A label at least6 is isolated in each required finite corner.
After one corner is omitted, the retained ring graph is a four-node path,
so at most two high contacts occur. Three low contacts determine the two
possible norm-two branches in the one-high case. Strict consecutive-cosine
bounds exclude every in-range branch. The two-high cases are excluded on
the entire enlarged real square by signed derivatives and a norm margin.
Auxiliary corner minors that enclose zero are not used for rejection.

This construction read is not a second implementation of the interval
engine. Its scope is the twenty literal rings and the specified topology;
it produces no new chamber. It is independent of the dimension29 theorem.
