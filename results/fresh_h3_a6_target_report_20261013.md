# Targeted H3+A6 completion and exact slice obstructions

This follow-up found new primitive reflective cuts in both recent H3+A6 hosts and built exact Coxeter cones with up to 773 walls. It did not find a finite-volume example. It also closes two specific low-dimensional search mechanisms without imposing a height bound: the original middle-weight slice contains no reflective normals, and an adjacent slice cannot remove a certified spacelike escape using any old-compatible real walls.

The starting lattices are Lα=<-α>+H3+A6 over O=Z[φ], α=φ or 2φ. Their complete height-30 prefixes have 45 and 56 walls. Their construction and limited determinant comparison are in the [previous report](fresh_h3_blocks_h9_report_20261013.md). No claim about worldwide novelty or nonreflectivity of the full hosts is made here.

## First assess the real geometry and integral lattice

Number fundamental weights starting at 1. Let

```
h = 2ω1(H3),
a = ω3(A6)+ω4(A6) = (1,2,3,3,2,1),
b = 7ω3(A6) = (4,8,12,9,6,3),
k = ω2(H3).
```

The saved A6-symmetry-fixed escapes lie in span(time,h,a). They have positive norm and satisfy every old wall inequality. Increasing their time coordinate by 1/100 gives an exact spacelike normal that cuts the original escape and has negative pairing with every positive-height old wall. Thus the continuous real cut region is present; the obstruction below is arithmetic.

The full integral intersection with span(time,ω1(H3),ω3(A6),ω4(A6)) has basis (time,h,a,b). This is not an assumed sublattice: the coordinate minor at indices (0,3,4,9) has an integral inverse. Its spatial norm is

```
6φ²s² + 6u² + 42uz + 84z²,
```

which lies in 3O for every integral s,u,z. The prime 3 is inert, so O/(3)=F9. The complete primitive candidate norms are divisors of 28, normalized by square units; none vanishes modulo 3. The root equation would require

```
-αx² = q  modulo 3.
```

But the norm of α from F9 to F3 is 2 for both time coefficients, whereas the norm of a nonzero rational q is 1. Hence -q/α is a nonsquare in F9. There is no solution.

Therefore **this entire four-dimensional field slice contains no primitive reflection normal**, including asymmetric combinations of the two middle A6 weights. A bounded symmetric-slice check with s,u≤1000 found no integral norm solutions before this stronger residue proof was archived. Raising that bound would add no information.

## A neighboring H3 direction produces actual cuts

Adding k=ω2(H3) breaks the modulo-3 obstruction. The new slice is

```
r = x·time + s h + t k + u a,
q(r) = -αx² + 6φ²s² + 4φ³st + 2φ²t² + 6u².
```

The integral basis (time,h,k,a) is also saturated, witnessed by the integral inverse of its coordinate minor at indices (0,1,3,4). Its dominant parameters satisfy s,u≥0 and t≥0. Reflection integrality excludes primitive norms divisible by 7 in this slice: it would force all four integral parameters to be divisible by 7. The remaining candidate pools are {1,2,4} for α=φ and {2,4} for α=2φ. In the former case norm 4 produces only imprimitive normals and is discarded.

One bounded search covered 0≤s,u≤100 and 0<t≤100. The t=0 slice is already covered by the residue obstruction. Exact conjugate bounds are

```
s'² ≤ q/2,   t'² ≤ 3q/2,   u'² ≤ q/6.
```

They follow from the diagonal entries of the inverse restricted positive Gram matrix. The search tested the exact norm equation, primitive integral reflection coefficients, all original wall inequalities, every Coxeter angle, and the common strict timelike center. It saved all compatible roots in this box, including those that do not cut the original target.

| α | New compatible roots | Roots cutting the original symmetric escape | Runtime |
|---|---:|---:|---:|
| φ | 772 | 734 | 12.26 s |
| 2φ | 362 | 237 | 13.86 s |

For α=φ, a particularly simple norm-one cut has parameters

```
(x,s,t,u) = (2+3φ, 0, 1+φ, 1+φ).
```

Its original-basepoint height is 21+34φ and its pairing with the old symmetric escape is (7-3φ)/19>0. For α=2φ, one norm-two cut has parameters

```
(x,s,t,u) = (3+4φ, φ, 2φ, 1+2φ),
```

height 40+65φ and escape pairing (4φ-2)/5>0. These are genuine primitive lattice reflections, individually compatible with every old wall. They are not finite-volume polytopes by themselves.

The full dictionaries are [α=φ](fresh_h3_a6_target_phi_neighbor_20261013.json) and [α=2φ](fresh_h3_a6_target_2phi_neighbor_20261013.json). This parameter box includes high original-basepoint heights, but it is a targeted slice search, not a uniform extension of the previous Vinberg cutoff.

## Compatible Coxeter cones and surviving volume escapes

All cross-pairings between dictionary roots were checked exactly. Acute pairs with forbidden finite angles are conflicts, just as positive pairings are. In particular the golden-field angle test distinguishes π/5 from the inadmissible angle 2π/5.

The dictionaries have 365 and 372 conflict edges. Five greedy orders produced five distinct maximal compatible selections in each host. “Maximal” means no additional dictionary root can be added; no maximum-cardinality claim is made.

| α | Wall counts of the five cones, including original walls |
|---|---|
| φ | 474, 773, 757, 746, 756 |
| 2φ | 153, 362, 347, 350, 331 |

Every selected cone has a strict timelike interior, valid Coxeter angles for every wall pair, and a certified spacelike extreme ray. Their artifacts are [α=φ](fresh_h3_a6_target_phi_cones_20261013.json) and [α=2φ](fresh_h3_a6_target_2phi_cones_20261013.json).

There is also an exact escape satisfying the entire finite dictionary of inequalities in each host. Their norms are 64/35-φ and (43882-26036φ)/6727, respectively. The union of all dictionary walls need not be Coxeter; these common escapes are used only to certify that every compatible subset of this finite dictionary still has infinite volume.

## The neighboring slice cannot close the volume at any bound

For the neighboring slice S=span(time,h,k,a), define the restriction of pairings by

```
P(w) = (-αw0, 2w1, w2, w6+w7).
```

The wall indices below are zero-based indices in the original 45- or 56-wall prefix. For α=φ, the common dictionary escape v has norm 64/35-φ>0 and satisfies

```
P(v) = (φ/2)P(w0) + P(w1)
       + (φ-3/2)P(w21) + (5/2-3φ/2)P(w33).
```

Every coefficient is strictly positive. For α=2φ, a saved escape v from the third selected cone has norm (62-38φ)/7>0 and satisfies the shorter identity

```
P(v) = (4-2φ)P(w5) + (φ-1)P(w12).
```

Again both coefficients are positive. For any real normal r in S that is acute-compatible with every retained original wall, either identity gives (r,v)≤0. The same spacelike escape therefore survives every such added wall, independent of its norm, arithmetic field, coefficient size, or search height.

Consequently **no full-dimensional finite-volume completion using only old-compatible real normals from this neighboring slice is possible**. This is stronger than a bounded failed search. It does not exclude cuts outside the slice, removal of original walls, or finite-volume polytopes elsewhere in either full lattice.

The [exact obstruction artifact](fresh_h3_a6_target_obstructions_20261013.json) includes both saturated bases, the modulo-3 census, the geometric cutters, both surviving escapes, and all positive projection coefficients.

## Independent verification and reproduction

The [independent audit](fresh_h3_a6_target_independent_audit_20261013.json) uses standard-library Fraction and integer arithmetic. It reconstructs the hosts and slice bases, checks all accepted primitive roots and reflection coefficients, every original/new and new/new Coxeter angle, the complete conflict graphs, maximality of all ten selections, strict interiors, every escape, and both global obstructions. It independently counts the rectangular embedding alphabets. The complete positive-height candidate search is source-reviewed; a second algorithm does not repeat every candidate triple.

Four focused regressions passed: the norm polynomial is compared with direct full-host pairings, both exact obstruction certificates are replayed, and forged accepted-root norms are rejected in each host.

```sh
PYTHONPATH=src:. python3 scripts/fresh_h3_a6_target_search.py --alpha phi --bound 1000 --output results/fresh_h3_a6_target_phi_search_20261013.json
PYTHONPATH=src:. python3 scripts/fresh_h3_a6_target_search.py --alpha 2phi --bound 1000 --output results/fresh_h3_a6_target_2phi_search_20261013.json
PYTHONPATH=src:. python3 scripts/fresh_h3_a6_target_neighbor.py --alpha phi --bound 100 --output results/fresh_h3_a6_target_phi_neighbor_20261013.json
PYTHONPATH=src:. python3 scripts/fresh_h3_a6_target_neighbor.py --alpha 2phi --bound 100 --output results/fresh_h3_a6_target_2phi_neighbor_20261013.json
PYTHONPATH=src:. python3 scripts/fresh_h3_a6_target_cones.py --alpha phi
PYTHONPATH=src:. python3 scripts/fresh_h3_a6_target_cones.py --alpha 2phi
PYTHONPATH=. python3 scripts/fresh_h3_a6_target_verify.py
PYTHONPATH=src:. python3 -m pytest -q tests/test_fresh_h3_a6_target.py
```

The two stated slices are now closed for their respective mechanisms. No further bound increase is warranted there. The full arithmetic hosts remain unresolved.
