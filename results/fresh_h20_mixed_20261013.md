# Mixed H20 constructions and finite-volume-preserving cuts

This round constructs a convex finite-volume H20 polytope from four
short and four long facets of Borcherds' H21 chamber. It has **767 true
facets**, including 415 projectively nonrational supporting planes over
`Q(sqrt(2))`. It is **not Coxeter**. The construction is a new branch of
this repository's search; no claim of mathematical priority is made.

The round then excludes two explicitly listed plane dictionaries, and
tests inward cuts whose finite volume follows from exact containment.
No finite-volume H20 Coxeter polytope is claimed. These finite searches
do not exclude other planes, other algebraic fields or other assemblies.

## The eight-tile mixed construction

The frozen input is the original 210-wall H21 chamber, with wall norms
2 and 4. The underlying example is classical; see the
[previous construction report](new_dimension_h20_facet_gluing_20261013.md)
for its primary sources and the distinction between short and long
sections. The new gluing uses their literally identical H19 intersections,
including the 81-facet and 50-facet seams that the earlier six-short-tile
patch did not use.

Let `a,b,c` be ambient walls **1,22,37**, with Gram

`[[2,-1,0],[-1,2,-2],[0,-2,4]]`.

These form a B3 diagram, with a short endpoint, a short central node,
and a long endpoint. Their intrinsic H20 ridge sectors are respectively

`gamma=acos(sqrt(2/3)), beta=pi/2-gamma, delta=pi/4`.

Let `M_ab` be reflection in `a+b`. For a short/long pair `s,l`, the
mixed transport `M_sl` is reflection in `l+sqrt(2)*s`. It fixes their
common seam pointwise and sends `l` to `-sqrt(2)*s`, hence maps the long
section into the short section's hyperplane on the opposite side of the
seam. Its coefficients belong to `Q(sqrt(2))`.

The eight source sections have parent labels

`a,b,c,c,b,a,c,c`.

Their consecutive attachment maps are

`M_ab, M_bc, reflection(a), M_bc, M_ab, M_ac,
 reflection(2b+c), M_ac`.

Both self-reflections here preserve the long section's hyperplane.
Their sector sequence is

`gamma,beta,delta,delta,beta,gamma,delta,delta`,

with total angle `2pi`. The seam facet counts around this cycle are

`208,50,81,50,208,81,50,81`.

The full transported wall sets verify all eight facet identifications,
and the last attachment returns the first complete tile exactly. Six
additional pairs of opposite normals do **not** give equal full facets
and are not counted as gluings. The separate two-tile mixed seam tests
give 264 and 263 uncancelled planes for the B2 and orthogonal interfaces.

The construction and all exact wall lists are in
[the script](../scripts/fresh_h20_mixed_b3_20261013.py) and
[the geometric archive](../data/enumerations/fresh_h20_mixed_b3_20261013.json).
Each stored normal is a pair of integer vectors `(A,B)`, denoting
`A+sqrt(2)B`; its sign is preserved under normalization.

## Global geometry is certified; the boundary is not Coxeter

The [support checker](../scripts/fresh_h20_mixed_support_20261013.py)
and [certificate](../data/enumerations/fresh_h20_mixed_support_20261013.json)
provide 6,136 exact nonnegative cone expressions: every one of the 767
boundary normals supports each of the eight tiles. Every coefficient
belongs to `Q(sqrt(2))`, with its sign checked exactly. Numerical LP only
selects possible cone generators. Pulling back to the rational source
section reduces coefficient recovery to exact rational linear systems.

There are also exact separating-plane certificates for all 28 tile
pairs. Thus the tile interiors are disjoint. Together with the eight
complete facet gluings and globally supporting remaining boundary
planes, these certify that the union is convex. Each remaining plane
contains a complete facet of one of its owner tiles, so all 767 planes
are true facets. Finite volume follows from the eight finite-volume
source tiles.

An explicit genuine ridge still has Gram `[[6,-2],[-2,6]]` and squared
cosine `1/9`. It comes from ambient walls 24 and 49 of tile A, and has
indices **24 and 96** in the new zero-based boundary ordering. The
certificate gives an integer point of square **-76** lying on exactly
these two boundary planes, with every other boundary inequality strict.
This is an actual non-Coxeter ridge, not merely a bad pairing of
nonincident supporting hyperplanes.

## Selection from new supporting planes

The first dictionary combines:

- the earlier 580 rational boundary planes;
- all 767 mixed boundary planes;
- the targeted mirror `W=3r_0+(2-sqrt(2))r_74`, where subscripts refer
  to the old 580 ordering.

After exact projective deduplication it contains **1,142 planes**, of
which 416 are nonrational, with 20,350 complete Gram conflicts.
This is genuinely outside the earlier closed 580-plane dictionary.
For example, new mixed wall 343 cuts the old 363-wall compatible
chamber's archived spacelike escape, although it conflicts with some
walls of that chamber.

The allowed elliptic squared cosines over `Q(sqrt(2))` are
`0,1/4,1/2,3/4,(2+sqrt(2))/4`, corresponding to labels 2,3,4,6,8.
A positive outward-normal pairing also conflicts. These are conditions
on the **complete Gram**, not a count of actual ridges. A Coxeter
intersection represented using some redundant inequalities still has
a compatible minimal facet set, which is itself a subset of the
dictionary; the selection therefore does not assume every listed plane
is a necessary facet.

Each exact spacelike vector `v` supplies a necessary clause

`at least one chosen wall r has (r,v)>0`.

Otherwise `v` satisfies all chosen inequalities, and perturbing it
toward the shared strict timelike center gives a strict spacelike
escape, proving infinite volume. The search reuses the old 580 single-wall
deletion certificates and the 18 compatible-chamber certificates,
recomputing every pairing against the new dictionary. They yield 589
distinct escape clauses.

For the 1,142-plane dictionary, just two clauses already contradict
compatibility: the old deletion escapes for indices 8 and 101 can be
cut only by their corresponding original walls. Those two walls have
Gram `[[2,-2],[-2,6]]`, with forbidden squared cosine `1/3`.
Thus **no subset of these 1,142 planes is a finite-volume Coxeter
polytope**. The exact two-variable UNSAT proof, vectors and all plane
coordinates are in the [dictionary archive](../data/enumerations/fresh_h20_mixed_dictionary_20261013.json);
[the script](../scripts/fresh_h20_mixed_dictionary_20261013.py) independently
replays the necessary clauses and Gram conflicts.

## Adaptive new mirrors: a bounded 1,262-plane dictionary

The [adaptive script](../scripts/fresh_h20_mixed_adaptive_20261013.py)
uses successive UNSAT cores to generate new planes, rather than raising
a previous search cutoff. For a rational blocking pair it rotates one
wall to an allowed angle with the other whenever the required coefficient
lies in `Q(sqrt(2))`, the common timelike center remains strict, and the
new plane cuts at least one current core escape. It considers both
orientations and labels 2,3,4,6,8.

The preset limit adds **120 distinct new mirrors**: 49 are positive
combinations supporting the original intersection, while 71 have an
inward coefficient and can trim it. The final dictionary has **1,262
planes**. Its final core contains 30 exact escape clauses and 6,699
conflict edges; replay confirms UNSAT. No compatible candidate reaches
the subsequent volume gate in this bounded round.

The [adaptive archive](../data/enumerations/fresh_h20_mixed_adaptive_20261013.json)
records every parent pair, exact coefficient, new plane, escape clause
and final core. This excludes all subsets of the **literal final
dictionary**, including its inward mirrors. It is not an exclusion of
inward cutting generally. The following sequential cuts generate
additional planes and impose exact containment explicitly.

## Inward cuts preserve finite volume automatically

Replacing a wall `r` by `v=r+t*s`, with `t<0`, while retaining `s`, gives

`r=v+(-t)*s`.

The two coefficients are positive. Thus the new inequalities imply the
old ones, and the new polytope lies inside the previously certified
finite-volume polytope. Keeping the common strict timelike point proves
nonemptiness and full dimension. This removes dependence on failing to
find spacelike escapes as a volume test.

For the old pair `r_0,r_74`, the first cut is

`v=3r_0-sqrt(2)*r_74`.

It has norm `48+24sqrt(2)`, makes angle `pi/8` with `r_74`, and pairs
with the strict center as `-48+10sqrt(2)<0`. The identity
`r_0=(v+sqrt(2)*r_74)/3` proves containment exactly.

The [inward script](../scripts/fresh_h20_mixed_inward_20261013.py)
checks redundancy after each replacement. It removes an inequality
only after recovering and checking an exact nonnegative combination
of the remaining field-valued normals. The saved replay anchors its
first 580 normals and timelike center to the authoritative source,
checks every new normal's positive norm, and checks all containment
and pruning identities.

The [12-cut archive](../data/enumerations/fresh_h20_mixed_inward_20261013.json)
contains twelve nonempty finite-volume regions. The first cut makes no
other inequality redundant and increases complete Gram conflicts from
7,072 to 7,485. A degree-based greedy choice then repeatedly projects
the new wall orthogonally toward walls 2 and 3. It removes no additional
walls and stalls at 7,536 conflicts, retaining 580 inequalities.

These counts are **raw pair conflicts**, not actual ridge counts.
The absence of a pruning certificate is also not used as proof that
every retained inequality is a true facet. The replay proves containment
and the reported exact Gram calculations.

There is a separate exact geometric failure certificate after **all
twelve cuts**. The original ridge at old boundary indices 23 and 100
survives: an integer point of square `-762753024` lies on exactly these
two retained inequalities and satisfies every other final inequality
strictly. Their Gram is `[[6,-2],[-2,6]]`, with squared cosine `1/9`.
Thus the final inward-cut polytope is **definitely non-Coxeter**,
independently of any unverified redundancy or pair-incidence assumption.
The alternating sequence is archived as a failed greedy choice. It is
not continued with a larger iteration cutoff.

## A block projection changes the required number field

Write `a=r_0`, `b=r_74`, and let `u=r_2`, `w=r_3`, again in the old
580 ordering. Their Gram is

`[[4,-4,0,0],[-4,6,-2,0],[0,-2,2,-1],[0,0,-1,2]]`.

Its determinant is -8. The **two-dimensional** block spanned by `u,w`
is positive definite; the whole four-root span is not.
The simultaneous projection is

`b'=b+(4/3)u+(2/3)w`,

orthogonal to `u,w`, with `(b')^2=10/3` and `(a,b')=-4`.
Instead of infinitely alternating two projections, one may make the
single inward cut `3a-sqrt(2)b'`. It has positive norm
`128/3+24sqrt(2)` and center pairing `-48+14sqrt(2)<0`.

This fixes the two right-angle conditions simultaneously but changes
the angle with b. For the inward family `v=a+t*b'`, `t<0`, its squared
cosine with b varies between `5/9` and `2/3`. The only Coxeter label in
that interval is **5**, which is unavailable over `Q(sqrt(2))`.

The exact label-five candidate has t equal to the negative real root
of

`1475t^4-7080t^3+10476t^2-4752t-324=0`,

approximately `-0.0599372348773694`. It obeys

`(-4+10t/3)^2 = 6*((3+sqrt(5))/8)*(4-8t+10t^2/3)`.

Its center pairing is `-16-14t<0`, and

`a=v+(-t)b+(-4t/3)u+(-2t/3)w`

proves finite-volume containment while retaining b,u,w. Thus this is a
concrete candidate in a different algebraic field, with three prescribed
angles fixed at once. The
[independent block report](fresh_h20_cusp_repair_inward_20261013.md),
[script](../scripts/fresh_h20_cusp_repair_golden_block_20261013.py) and
[exact archive](../data/enumerations/fresh_h20_cusp_repair_golden_block_20261013.json)
use the positive scalar multiple `V=3a+t*(3b+4u+2w)`. They isolate t
within the rational interval `(-3/50,-59/1000)` and refine it by exact
bisection. This V has norm `36-72t+30t^2` and center pairing `-48-42t`.

Against the 579 retained original inequalities, its declared finite-label
screen finds 361 compatible pairs, 18 positive pairings and 200 other
intersecting angles. These are pair tests, not assertions about true
facets or actual ridges; no pruning was performed on this candidate.

More decisively, the same genuine old ridge 23/100 survives the quartic
cut. Its strict timelike ridge point has new-wall pairing
`-209088*(1+t)<0`. The two original bounding walls still have Gram
`[[6,-2],[-2,6]]`. The quartic polytope is therefore a verified nonempty
finite-volume **non-Coxeter** polytope, despite its three repaired local
angles. The unchanged ridge gives a concrete next repair target; it
does not obstruct other simultaneous block cuts.

## Replays

```sh
python3 scripts/fresh_h20_mixed_b3_20261013.py
python3 scripts/fresh_h20_mixed_support_20261013.py --verify
python3 scripts/fresh_h20_mixed_dictionary_20261013.py --verify
python3 scripts/fresh_h20_mixed_adaptive_20261013.py --verify
python3 scripts/fresh_h20_mixed_inward_20261013.py --verify
python3 scripts/fresh_h20_cusp_repair_golden_block_20261013.py --verify
```

The first command rebuilds exact tile transports and seam identifications;
the others replay archived certificates without repeating LP searches.
All five mixed-construction scripts pass compilation. Exact support,
dictionary, adaptive-mirror, inward-containment and quartic-block replays
pass. The mathematical endpoint remains a
finite-volume non-Coxeter mixed patch and explicitly scoped construction
obstructions, together with concrete inward cuts for further work.
