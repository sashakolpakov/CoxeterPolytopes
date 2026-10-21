# Single-factor cusp sections: leaf sums and all rank-two factors

No Coxeter H20 section of the archived Borcherds H21 chamber has its normal in any of the 79 leaf-sum planes constructed here. The exclusion covers every real shift in each plane. The same calculation closes the entire real normal span of each of the four rank-two affine factors in the complete 18-cusp census: two factors of type C2 and two of type A2. These statements overlap in two leaf-sum planes; they are not counts of disjoint search families.

The 79 planes are distinct from each other and from all 59 previously excluded twin-leaf **difference** planes. Their local cusp links are Coxeter throughout the relevant open parameter interval, so the failures occur at certified actual ridges elsewhere in the section. The result leaves arbitrary directions in the higher-rank single factors open. It does not exclude all H20 sections or H20 Coxeter polytopes.

The [standard-library Fraction checker](../scripts/fresh_h20_single_factor_links_20261013.py) and [exact certificate](fresh_h20_single_factor_links_20261013.json) reconstruct the original 210 wall normals, all 18 cusp rays and all 53 positive affine relations. They check 156 exact ambient face points. No floating-point predicate, coefficient bound, or finite angle-label cutoff is used.

## Actual-face certificates

Write the original Lorentz form as `Q`, the outward wall normals as `r_i`, and the saved strict center as `h`, with `h²=-14`. Every face used below has a positive-definite wall Gram matrix `A` whose inverse is entrywise nonnegative. The checker projects `h` orthogonally to those walls and verifies the resulting point is future timelike, has exactly the declared zero wall pairings, and is strictly negative against every other one of the 210 walls.

Positive combinations of these future timelike points remain future timelike. The zero set of such a combination is exactly the intersection of the endpoint zero sets. The interpolation arguments below therefore prove actual section ridges, rather than inferring adjacency from a projected two-wall Gram matrix. The two projected ridge normals are independent in every case, as shown below; a small perturbation within the section makes both inequalities strict. Thus these are genuine interior sections.

## Seventy-nine leaf-sum planes

In one affine factor choose two or three mutually orthogonal norm-two leaves with a common neighbor and identical pairings outside the chosen leaves. Their positive affine marks are equal. Let `S` be this set, `k=|S|`, and normalize the positive null relation to

`Delta = sum_{i in S} r_i + sum_{j outside S} m_j r_j`, with every `m_j>0`.

Set

`e0 = sum_{i in S} r_i`,  `e(t)=e0+t Delta`,  `e(t)²=2k`.

The literal census gives 59 two-leaf choices and 20 three-leaf choices. The latter come from the five affine D4 factors at cusp 58. All other cusp factors are orthogonal to `e(t)` and remain unchanged.

Projection of the selected roots gives the positive affine relation `sum r_i'=-t Delta`. The complementary roots give `sum m_j r_j'=(1+t)Delta`. Hence the full local link is bounded and Coxeter precisely for `-1<t<0`. In 77 cases the two groups are orthogonal affine simplex factors, retaining every original factor inequality. The checker verifies their Gram matrices, positive kernels, positive-definite proper minors, Coxeter pairings, and total quotient rank `n-1`.

In the two rank-two C2 cases the complementary root instead projects to `(1+t)Delta`. It is a strict redundant inequality on the horosphere, not a facet. The selected pair is an affine A1 factor. For `t>=0` or `t<=-1`, the normal is a one-sign combination of the old roots and its hyperplane cannot meet the strict chamber interior.

Every one of the 79 constructions has a norm-two auxiliary wall `j` attached with pairing `-1` to exactly one selected leaf and orthogonal to every other root of the affine factor. Choose a different selected leaf `i`. Then

`(e0,j)=(Delta,j)=-1`,  `(e0,i)=2`,  `(Delta,i)=0`,  `(i,j)=0`.

The projected pairing is consequently

`(i',j')=(1+t)/k > 0`.

To prove it is an actual obtuse ridge, the checker supplies two exact face points. The left face has all complementary factor roots and `i,j` zero; the right face has all selected leaves and `j` zero. Their zero sets intersect exactly in `i,j`. After positive rescalings these points, called `X,Y`, satisfy

`(e0,X)=(Delta,X)=(Delta,Y)=-c`,  `(e0,Y)=0`, with `c>0`.

Thus `x(t)=(-t)X+(1+t)Y` lies in `e(t)`-perp for every `-1<t<0` and has exactly the two required zero inequalities. The projected ridge Gram determinant, with `u=1+t`, is

`4-4/k-u²/k > 0`.

This certifies a genuine obtuse ridge throughout the entire interval. No Coxeter label can repair it.

For distinctness, each normal plane `span(e0,Delta)` has the unique null line `span(Delta)`. Planes with different null lines differ. At the same cusp, adding a multiple of `Delta` does not change the normal's pairings with incident walls. The new sum planes have three or four nonzero incident pairings, whereas every archived difference plane has two. Among the sum planes those supports differ. The checker verifies these tests literally, including positivity and null-orthogonality of the old plane generators.

## Complete C2 factor spans

The two factors are `[2,27,184]` at cusp 1 and `[0,32,162]` at cusp 2; all wall indices in this report are zero-based. In each factor label the norm-two roots `a,b` and the norm-four root `c`. Their Gram matrix is

```text
 2  0 -2
 0  2 -2
-2 -2  4
```

With `Delta=a+b+c`, every spacelike normal in the factor's full three-dimensional span is uniquely

`e=x a+y b+t Delta`,  `(x,y)!=(0,0)`,  `e²=2(x²+y²)`.

An interior section is possible exactly when

`-max(x,y,0) < t < -min(x,y,0)`.

Necessity follows by normalizing a strict point to `(Delta,p)=-1`: the negative factor pairings are three positive numbers summing to one. The constructions below also establish sufficiency. Equality at either endpoint cannot meet the strict interior.

For each factor the certificate finds orthogonal short roots `j,k` such that `j` pairs only with `a`, with value `-1`, and `k` pairs only with `b`, also with value `-1`. The three ambient faces formed by `j,k` and two of `a,b,c` are positive definite. Normalize their exact points to `(Delta,p)=-1`. Their factor pairings are respectively `(-1,0,0)`, `(0,-1,0)`, and `(0,0,-1)`.

For every strict parameter interval, choose strictly positive barycentric weights with weighted mean of `x,y,0` equal to `-t`. Their combination is a timelike section point whose only zero inequalities are `j,k`. These projected normals are independent: otherwise `e` would be in their span and would vanish on all three face points, forcing `x=y=t=0`. Thus they form an actual ridge with projected pairing

`(j',k')=-(x+t)(y+t)/(2(x²+y²))`.

If `xy<=0`, the two factors in the numerator have opposite signs throughout the strict interior interval. The actual ridge is obtuse. This includes `x=0` or `y=0`.

If `xy>0`, reverse the overall orientation and swap the labels if needed so that `x>=y>0`. Now `-x<t<0`. Use the exact faces `{b,c,j}` and `{a,b,j}`, normalized as above. The positive combination with coefficients `-t` and `x+t` is a timelike section point with exactly `b,j` zero. Independence follows because `e` has nonzero values on both endpoint face points. Its projected ridge pairing is

`(b',j')=2y(x+t)/(2(x²+y²)) > 0`.

The checker also supplies the swapped case. This covers all real spacelike normals in both factor spans, not merely the leaf-sum planes or a finite list of local angles. A pure multiple of `Delta` is null and does not define a hyperbolic section meeting the strict interior.

## Complete A2 factor spans

The two factors are `[152,177,207]` at cusp 7 and `[152,177,209]` at cusp 390. Every root has norm four and every distinct pair has pairing `-2`. A normal whose hyperplane meets the strict interior must have both coefficient signs in this affine root basis. Reversing the normal if necessary, write it as

`e=x r_i+y r_j-t r_k`,  `x,y>=0`,  `x+y>0`,  `t>0`.

This includes the zero-coefficient boundary cases. Its square is

`q=4(x²+y²+t²-xy+xt+yt)`,

and

`q-4t²=4((x-y)²+xy+t(x+y))>0`.

For each of the three choices of negative root in each factor, the checker finds short walls `a,b` with `a.b=-1`, where `a` is orthogonal to the entire factor, `b` is orthogonal to the two positive roots, and `b.r_k=-2`. The exact ambient faces `{a,b,r_i,r_j}` and `{a,b,r_k}` are positive definite. The normal has opposite strict signs at their certified face points. A positive interpolation therefore lies in the section, is timelike, and has exactly `a,b` zero. Independence again follows from the nonzero endpoint evaluations.

The squared cosine at this actual ridge is

`q/(4q-8t²)`, strictly between `1/4` and `1/2`.

There is no integer Coxeter label in that interval. All six negative-root cases are independently reconstructed. The earlier affine-A2 argument covered an orbit represented by cusp 390; the present source also checks the other factor at cusp 7 directly. No symmetry transfer is needed for this certificate.

## Replay and remaining scope

Run:

```sh
python scripts/fresh_h20_single_factor_links_20261013.py --verify
```

The continuous arguments are the proofs above; the checker certifies their exact Gram, affine-relation, face, interpolation, and distinctness data. It is not a formal proof assistant. In particular, the rank-two all-parameter claims do not come from sampling, and the actual-ridge claims do not come from hyperplane intersection tests alone.

The final saved-result replay passed. A [separate SymPy reconstruction and proof review](fresh_h20_single_factor_links_20261013_review.md) checks all 156 ambient face points against all 210 walls, the 79 symbolic ridge identities, and every rank-two auxiliary configuration without importing this producer's arithmetic. The continuous coverage arguments are reviewed separately from that exact reconstruction.

Together with the earlier mixed-factor and pure-A1 results, this narrows the remaining cusp-section search to higher-rank single-factor directions outside the excluded planes. It does not classify all such directions, sections avoiding the listed cusp orbits, changed ambient chambers, or reflection subdivisions of these newly studied sections.
