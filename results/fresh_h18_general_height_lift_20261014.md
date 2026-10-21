# Arbitrary real heights do not complete the bare H18 lift

Every lift in the class below has infinite volume in every dimension at least 20. This removes the two-axis, equal-length and symmetric-height restrictions of the earlier H18-to-H20 experiment. It does **not** exclude adding walls with different source projections, changing the retained anchors, or other dimension-raising constructions.

The [exact certificate](fresh_h18_general_height_lift_20261014.json) and [checker](../scripts/fresh_h18_general_height_lift_20261014.py) reconstruct the source Gram inertia, common center, six fixed anchor contacts and three null-ray sign lists with independent Fraction arithmetic over Q(sqrt(2)). Generation and saved replay pass. The [separate arithmetic and proof review](fresh_h18_general_height_lift_20261014_review.md) also passes, independently checking the six anchors and all 160 center/null-wall signs. The universal argument uses these literal identities; it does not sample heights or truncate the angle labels.

## The broader class being tested

Let the forty source normals be those of the [archived Vinberg H18 hybrid](../data/enumerations/vinberg_h18_mixture_facets_2026_10_13.json). Its Gram has inertia `(18,1,21)`, so the nondegenerate source space W has signature `(18,1)`. There are 26 norm-two roots and fourteen norm-one roots. The sum h of all forty roots is timelike and pairs strictly negatively with every root, as checked directly from the full Gram.

Adjoin a positive Euclidean space U of dimension k, with k at least two. Consider a finite collection of proposed Coxeter walls in W orthogonal-sum U satisfying these conditions:

- Retain all 26 norm-two source normals with zero U-component.
- Every other normal, after positive rescaling, has the form `r_i+v`, where r_i is one of the fourteen norm-one source roots and v is any real vector in U. Multiple heights over the same root are allowed.
- Over each of the three roots `L18`, `L25`, `L26`, the height set is nonempty and its convex hull contains zero. This includes asymmetric sets and is the condition that convexly averaging these lifted normals recovers the source normal.
- Every pair of distinct proposed walls satisfies the acute Coxeter condition.

There is no norm bound, coefficient field, prescribed height direction or symmetry assumption. Even convex recovery over the other eleven short roots is unnecessary. The point `(h,0)` is automatically strict timelike for every wall.

The [old two-axis construction](fresh_dimension_construction_routes_20261014.md) is a special case before its six cusp-repair walls. Those later repairs have source projections proportional to sums of source roots and are outside the class above. Its already certified 60-wall chamber is therefore not being subjected to a contradictory new claim.

## Fixed anchors force very small core height sets

Each core root has two retained norm-two anchors with pairing minus one:

| Core root | Retained anchors |
|---|---|
| L18 | L17, L21 |
| L25 | L1, L9 |
| L26 | L5, L13 |

For a lift `r_i+v`, the squared cosine against either anchor is

```text
1 / [2(1+|v|^2)].
```

It is strictly positive and at most 1/2. Among Coxeter angles, only pi/3 and pi/4 have squared cosines in that interval: they give `|v|^2=1` and `|v|^2=0`, respectively. In particular, arbitrary large labels or dotted contacts cannot evade this anchor constraint.

For two distinct heights v,w over the same core root, acuteness requires `1+v·w <= 0`. Since both lengths are at most one, equality in Cauchy–Schwarz is forced: both are unit vectors and w=-v. There cannot be a third distinct height. Together with zero in the convex hull, each core height set is consequently either `{0}` or an antipodal unit pair. This deduction establishes symmetry rather than assuming it.

The three core roots are pairwise parallel in the source: `r_i·r_j=-1`. For two antipodal unit height pairs, put `c=v_i·v_j`. Their four cross-wall relations give normalized nonnegative weights `(1-c)/2` and `(1+c)/2`. Every nonzero finite Coxeter weight is at least 1/2. Hence either c=0 or |c|=1: the two height lines are orthogonal or equal. If one weight is 1, that relation is parallel; it is included.

## Some source cusp always has an unused added direction

When k=2, a core zero-height set immediately gives a pair of core roots whose combined heights span at most one dimension. Otherwise there are three unit lines, pairwise equal or orthogonal, in a two-plane. At least two lines are equal. Thus again some pair of core roots has combined height span less than k.

When k is at least three, each core contributes at most one line. Any pair has height span at most two, which is already less than k. Choose such a pair i,j and a unit vector u perpendicular to all its heights.

For each of the three possible pairs the certificate verifies, without floating point,

```text
delta = r_i+r_j;
delta^2 = 0,   delta·h < 0;
delta·r_i = delta·r_j = 0;
delta·r_a < 0 for every other norm-one source root;
delta·r_b <= 0 for every norm-two source root.
```

The strict inequalities for all other short roots are essential: no arbitrary finite heights over those roots can close the cusp transversely at delta.

## A strict ambient escape

Set `tau=-1/(4 delta·h)>0` and, for sufficiently small positive epsilon, take

```text
X = delta + epsilon u + epsilon^2 tau h.
```

Orthogonality of W and U gives exactly

```text
X^2 = epsilon^2/2 + epsilon^4 tau^2 h^2 > 0,
X·h = delta·h + epsilon^2 tau h^2 < 0.
```

The two selected core families pair strictly negatively with X: their delta and u terms vanish, and their h term is negative. Every retained long wall has nonpositive delta pairing, no u term and a negative h term. Every other short-root lift has a strictly negative constant delta term; its arbitrary finite height contributes only a term of order epsilon. Since the collection is finite, one sufficiently small epsilon makes all those inequalities strict simultaneously.

Thus X is a future-oriented strict spacelike point of the chamber cone. Together with the strict timelike point h, it yields an open patch at infinity and infinite hyperbolic volume. The proof does not infer ambient volume from a lower-dimensional face.

This excludes the entire stated class in H20 and in every H(18+k) for k at least two. Trying non-axis directions, unequal lengths, additional copies of the same short-root directions, or asymmetric finite height sets cannot fix the bare lift while retaining its hypotheses. The constructive alternatives must introduce new source projections, replace the fixed anchors or change the core recovery condition. The result does not supply a new finite-volume example or exclude those alternatives.

```sh
python3 scripts/fresh_h18_general_height_lift_20261014.py --verify
```
