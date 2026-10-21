# Independent review of all last-block 4/5 contact profiles

**PASS.** The [independent checker](../scripts/fresh_h9_star_lastblock_20261014_review.py) reconstructs every new domain restriction and all **77,952** residual determinants using standalone Fraction arithmetic in the basis `1,sqrt(2),sqrt(5),sqrt(10)`. All 32 zero determinants are excluded by independently repeated actual-vertex LDL witnesses. The [review certificate](fresh_h9_star_lastblock_20261014_review.json) binds the frozen [producer](../scripts/fresh_h9_star_lastblock_20261014.py) and [data](fresh_h9_star_lastblock_20261014.json). Independent generation and saved replay pass.

This excludes the specified case-0 branch: all six leaves of the first three blocks have own contacts 3/3, their internal finite labels are arbitrary, and the last internal label is **4 or 5**, with all four last own contacts arbitrary finite Coxeter labels. It does not address other last internal labels or the whole incidence.

## New all-label reduction

All six retained-own-root budget supports and all four triangles involving the third root 14, a last anchor, and one last-block root are checked against literal vertices of the 378-vertex incidence.

Suppose an anchor's own cosine alpha has `alpha²>=3/4`, corresponding to label at least 6. The retained-root budget forces every foreign contact right: its retained own contribution is at least 3/4 and a foreign nonright contact costs at least 1/2. In each literal triangle with root 14, either additional nonright bond would already give a nonpositive three-node Gram, since `3/4+1/4=1`. Hence the anchor is right to root 14, and root 14 is right to that high-contact support node.

Let c be the last internal cosine, `D=1/(1-c²)>=2`, beta the other nonnegative own cosine, `X=D(alpha²+beta²+2c alpha beta)`, and Q the projection cost of root 14 in its **own** first block, excluding any last-block contribution. Thus `Q>=1+sqrt(2)/2`. If root 14 has no last-block foreign contact, the two residual normals are orthogonal and both negative. If it contacts the other last-block node at pi/3, their determinant is

`(X-1)(Q-1) + D(alpha²-1)/4`.

The independent checker verifies the coefficient identity using `D(1-c²)=1`. Its lower bound is strictly positive. More explicitly, put `a=alpha²-3/4>=0` and `d=D-2>=0`; after minimizing beta and Q, the lower expression becomes

`sqrt(2)/4-1/8 + d(3sqrt(2)/8-1/16) + (2+d)a(sqrt(2)/2+1/4)`.

Every displayed coefficient is positive. Therefore the residual contains a negative definite plane, impossible in Lorentz signature (9,1). This removes every own label at least 6 without a cutoff or a numerical inequality.

## Finite domain and literal witnesses

The first three internal labels are forced at least 4 by their missing triples with own contacts 3/3. Three explicitly checked placements of the previously reviewed two-dense theorem exclude labels at least 6. These inherited theorem proofs are not repeated here. All 24 cross-block right-angle guards are checked using their literal containing vertices and the positive null vector of the relaxed (4,3,4) path.

The remaining last-block own profiles are independently enumerated by the exact hyperbolic-triangle reciprocal condition: 11 ordered profiles when the internal label is 4 and 13 when it is 5. All 290 ordered profile pairs are covered. They split into 194 nonpositive own-triangle cases, ten inherited opposite-single cases, two inherited all-own-three cases, and 84 active cases (36 and 48 respectively). The inherited exclusions are bound through their frozen independent reviews; no duplicate replay of those old large censuses is claimed.

For each active pair, a common own root forces the two anchors right via a literal three-node all-ones witness. Both orientations of each retained-root budget are checked. A profile other than 3/3 has no foreign contacts; a 3/3 profile has at most one, necessarily pi/3. The checker independently reconstructs all 168 new actual-vertex label-4 failures that bound the two third-anchor contacts by 3. Each witness contains the tested pair. Perron monotonicity justifies relaxing every unlisted meeting contact to right in these failures.

## Complete arithmetic replay

The checker reconstructs all 16 positive eight-root block inverses, every allowed foreign pattern, and all residual matrices directly from their contact vectors. It computes each 3-by-3 determinant using the six signed permutations, rather than the producer's expanded determinant formula. Every saved coefficient is compared exactly.

All 19,488 anchor determinant occurrences are negative. Of the 77,952 complete three-root cases, 77,920 have nonzero determinant. Their total Gram would have rank at least 11, incompatible with H9. The remaining 32 cases have a nonpositive Gram at an actual vertex; every saved LDL pivot and the first failure order is independently checked. The complete record Cartesian product and empty unresolved list are verified.

No producer arithmetic, symbolic algebra library, finite-diagram classifier, metric optimization, or floating-point tolerance is used by this review. Its inherited dependencies are the already reviewed geometric incidence, two-dense gate, all-own-three family and opposite-single family, each used only where its literal hypotheses are checked.
