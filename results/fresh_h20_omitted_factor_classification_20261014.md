# One omitted facet of a single affine cusp factor

No new finite-volume polytope was found. For every affine factor of rank at least four at the eighteen saved Borcherds H21 cusp representatives, this certificate classifies all real section directions that retain exactly one fewer local facet and rejects every resulting H20 Coxeter section. The classification includes 48 irrational directions over Q(sqrt(6)); these have valid local Euclidean Coxeter links but actual obtuse ridges elsewhere.

The [local certificate](fresh_h20_omitted_factor_classification_20261014.json) covers 45 factors, 387 choices of omitted facet, and 93,921 possible contact supports. There are 175 accepted omission/direction pairs, representing 119 spatial planes modulo the cusp ray. Of these, 112 pairs lie in 56 previously excluded twin-leaf planes. The [global certificate](fresh_h20_omitted_factor_classification_20261014_global.json) treats the remaining 63 planes: 48 irrational and 15 rational. It supplies 52 continuous actual obtuse-ridge proofs, ten continuous actual angle-gap proofs, and one family reduced to a single algebraic parameter and rejected there. A second isolated vertex certificate covers an endpoint of another family.

This result does not classify the three affine rank-three factors. Rank-one and rank-two cases have earlier separate proofs. It assumes one active affine factor and exactly one omitted local inequality; it does not address arbitrary sections avoiding the saved cusps or changed ambient chambers.

## Complete local reduction

Let an active affine factor have rank n and n+1 roots. Write its positive null relation as `D=sum m_i r_i`, where D is a positive multiple of the cusp ray. After deleting the proposed omitted root, the n retained roots have positive-definite Gram B. Their diagram is a crystallographic forest. Normalize the spatial section normal to square one and put `a_i=(e,r_i)` on the retained roots. The projected Gram is

`H=B-aa^T`.

The retained local section is bounded, has dimension n−1 and n facets, and hence is a Euclidean Coxeter simplex. Thus H is irreducible affine, has a positive kernel, and has rank n−1. Conversely, these conditions certify its local Gram. The rank condition is exactly `a^T B^-1 a=1`, and its kernel is `k=B^-1 a`, oriented so that every k_i is positive.

For n≥4, an irreducible affine diagram has no triangle and its finite edge labels are only 3 and 4. Contacts of the same sign produce a complete negative-pairing graph, so there are at most two of either sign. Every opposite-sign pair needs an old edge to avoid a positive projected pairing. Since the old graph is a forest, two contacts of each sign are impossible. Therefore a has at most three nonzero entries. A mixed three-contact support must be an old three-node path.

At least one supported old node is adjacent to a node outside the support. Otherwise the support is a union of old components; the nonempty unsupported remainder is disconnected from it even after projection, contradicting irreducibility. At such a boundary edge, projection strictly increases the squared cosine. An old label-four edge cannot survive. An old label-three edge must become label four, so its supported node i satisfies `a_i^2=B_ii/2`. Overall sign reversal allows the positive square root.

For every other supported node j and every target squared cosine `c` in `{0,1/4,1/2}`, its contact must solve

`[a_i^2+c(B_ii-a_i^2)] a_j^2 - 2 B_ij a_i a_j + B_ij^2-c(B_ii-a_i^2)B_jj=0`.

The discriminant is rational. For the literal norm-two/norm-four source, every solution lies in Q(sqrt(2),sqrt(3),sqrt(5)). The checker enumerates both real quadratic branches, including a zero discriminant, and then verifies the full projected Gram, the unit condition, and the strict positive kernel. This is a finite consequence of the angle equations, not a coefficient grid or a sampled-angle cutoff.

There are 14,707 different exact local systems and 56,320 quadratic branch assignments before repetition across literal supports. All 175 accepted directions have support one or two; no support-three direction survives. All field equalities use the eight independent square-root basis elements. Signs use exact rational square-root enclosures, increasing precision when necessary; floating-point approximations are not part of the certificate.

## The whole shift interval and the omitted inequality

Choose the representative `e0=sum k_i r_i`, with `e0^2=1`. Every normal in the resulting spatial direction is proportional to

`e(t)=e0+tD`.

On the horosphere normalized by `(D,x)=-1`, let `y_i=-(r_i,x)`. The retained simplex has `y_i>=0` and `sum k_i y_i=-t`. Its omitted inequality is `sum m_i y_i<=1`. Testing the simplex vertices gives the complete interval

`-L<=t<0`, where `L=min_i(k_i/m_i)`.

Outside this interval the omitted wall is not redundant or the section does not meet the interior. At the lower endpoint, coincidence or contact with an existing local face is permitted; it is not silently discarded. In ambient coordinates,

`e0-LD=sum (k_i-Lm_i)r_i-Lm_o r_o`.

All displayed retained coefficients are nonnegative. If all vanish, the endpoint hyperplane is the old omitted facet and does not meet the strict interior of the source. Otherwise this identity proves global redundancy of the omitted inequality on that endpoint section. For interior values of the interval all retained coefficients in the analogous identity are positive.

## Global actual ridges, not just projected pairings

Write `t=-uL`, so `0<u<1`. For each of the 63 new planes the certificate gives two rational future timelike ambient face points, normalized by `(D,x0)=(D,x1)=-1`, with

`(e0,x0)=0`, and `(e0,x1)=-L`.

They come from exact projection of the strict source center onto positive-definite ambient faces. Every one of the 210 source inequalities is checked, including the exact zero set. Thus

`x(u)=(1-u)x0+u x1`

is future timelike, belongs to the section, and has exactly the common zero set K for `0<u<1`. The old roots in K are independent and positive definite. Their section projections are independent too: `e(u).x0=uL` is nonzero, whereas every root in K is orthogonal to x0. Since the projections lie in the positive-definite tangent space at x(u), their Gram is positive definite. A small perturbation inside this local orthant leaves any prescribed two of these walls active and makes the others strict. The reported pairs therefore bound genuine section ridges.

At u=1, the only additional zero wall is the omitted root. The nonnegative identity above makes it redundant among the common retained walls. This establishes endpoint actuality whenever the section meets the source interior. The exact angle inequalities are checked there as well; the sole remaining endpoint case has a separate actual finite-vertex certificate.

The projected pairings and squared norms are polynomials of degree at most two in u. The checker proves each necessary strict sign on the entire open interval by exact Bernstein coefficients, removing and recording any endpoint zero factors first. In the angle-gap cases it checks both squared norms and the cross-multiplied numerator/denominator inequalities; no numerical angle approximation is used.

Nine gap cases lie between consecutive Coxeter values `0,1/4,1/2`. One case, row102, has squared cosine strictly between `1/2` and `3/4` and uses the separate [all-angle source-section theorem](fresh_h20_crystallographic_angle_gate_20261014.md), which excludes label5 and every other finite label outside `{2,3,4,6}` in a hypothetical finite-volume Coxeter section. Its source manifest is rechecked by this certificate. This is not an arithmetic assumption about the proposed section.

The remaining family is row110, at cusp14, with omitted wall168. Its actual ridge `{2,188}` has

`cos^2(theta)=u^2/(4-u^2)`.

For `0<u<1`, only the Coxeter value `1/4` is possible, forcing `u=2/sqrt(5)`. At that exact value, a separately reconstructed finite vertex has a genuine ridge `{12,3}` with squared cosine strictly between zero and `1/4`. The endpoint u=1 is an old-facet hyperplane. The other isolated certificate is row83 at u=1, whose actual ridge `{4,1}` has the same forbidden open angle interval. Both vertices are checked against all 210 original inequalities, with exactly twenty active walls and an independent projected Gram.

The source is finite volume, and interior hyperplane sections inherit finite volume: in the Klein model they are sections of a finite polytope with finite or ideal vertices. The obstruction here is Coxeter ridge geometry, not failure to prove boundedness of a proposed new chamber.

## Reproduction and limits

Both generators and both saved-result replays pass on the frozen source files.

The [local independent review](fresh_h20_omitted_factor_classification_20261014_local_review.md) checks the completeness proof and replays the producer. Its separate Fraction/Q(sqrt(6)) implementation verifies all 175 accepted records and their shift intervals; it does not repeat the exhaustive support census.

The [global independent review](fresh_h20_omitted_factor_classification_20261014_global_review.md) separately reconstructs all 126 rational face points, 63 interpolation families, 208 polynomial positivity certificates, endpoint redundancy identities, and both exceptional vertices. Its standard-library arithmetic imports no producer code. Both independent saved replays pass; the older plane exclusion and the all-angle theorem remain explicitly identified dependencies.

```sh
python scripts/fresh_h20_omitted_factor_classification_20261014.py --verify
python scripts/fresh_h20_omitted_factor_classification_20261014_global.py --verify
```

The [global recipes](fresh_h20_omitted_factor_classification_20261014_global_recipes.json) specify only witness wall indices and the two exceptional parameters. The verifier reconstructs the exact points and identities; the numerical searches that suggested those indices are not premises.

The completed result covers every rank≥4 single-factor section with exactly one omitted local facet at these cusp representatives. The independent rank-three task produces further algebraic directions and is deliberately outside this result. No conclusion about general H20 nonexistence or a new polytope is drawn.
