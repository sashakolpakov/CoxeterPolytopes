# New seed metrics and simultaneous cap families

This numerical search changes the starting H9 metric instead of adding more
walls to the previously obstructed split seed. It finds no finite-volume
candidate. These are bounded discovery calculations, not exact exclusions.

There are eight pairs of seed walls, partitioned into two groups. If Z is
the eight-by-two group-incidence matrix, set

\[
K=Z(aJ-bI_2)Z^{\mathsf T},\qquad
G=\begin{pmatrix}2I_8&-2I_8-K\\-2I_8-K&2I_8\end{pmatrix}.
\]

The tested range has a>b+2 and 0<b<4/max(n1,n2). Each constructed Gram
matrix has numerical rank10 and signature(9,1). The old equal-group
construction occurs in this family after a parameter change; unequal
groups5+3,6+2,7+1 give different symmetry types.

For each prescribed finite pairing profile, kernel consistency requires
the sum of the two pairings at a matched seed pair to be constant within
each group. The existing block formulas enumerate these profiles modulo
permutations. The norm-two equation gives a curve in(a,b).

The first scout fixes seven rational b values and retains admissible roots
of every profile equation, giving388 parameter configurations in the three
unequal partitions. Every whole cap bank has a numerical spacelike escape.

The next scout searches intersections of two norm curves, allowing several
cap families to occur at the same metric. It uses500 samples in b and
refines sign changes; it can miss tangencies, close crossings, and roots
outside the stated a<20 window.

| Pairing labels | Sampled parameter intersections | Whole banks with no numerical outside vertex | Maximal compatible selections checked | Smallest selected maximum Klein square |
|---|---:|---:|---:|---:|
|2,3,5|3,415|2|188|15.3709|
|2,3,4,5|41,934|54|47,195|5.8284|

Every checked compatible selection has a vertex outside the Klein sphere.
The four-label configurations permit at most six mutually compatible caps;
the three-label configurations permit at most two. These maxima concern
only the banks in the table's third column. The whole-bank checks include
mutually incompatible caps: their apparent boundedness is only a necessary
screen, not a Coxeter construction.

Whole-bank geometry uses a sector for simultaneous permutations inside
each group. All banks contain the full prescribed cap orbits. Every point
can be permuted into this sector without changing its Lorentz norm.
Selected sets need not preserve the symmetry, so their geometry is always
recomputed in the full nine-dimensional Klein chart. The compatibility
screen allows finite labels2 through120 and dotted pairs. Every geometry
and compatibility decision is numerical. One whole-bank calculation used
Qhull's perturbation fallback; it supplies no exact boundedness claim.

Reproduction:

```sh
python scripts/fresh_h9_asymmetric_blocks_20261013.py
python scripts/fresh_h9_asymmetric_intersections_20261013.py
python scripts/fresh_h9_asymmetric_intersections_20261013.py --select
python scripts/fresh_h9_asymmetric_intersections_20261013.py --labels 2 3 4 5 --output data/enumerations/fresh_h9_asymmetric_four_labels_20261013_scout.json
python scripts/fresh_h9_asymmetric_intersections_20261013.py --select --output data/enumerations/fresh_h9_asymmetric_four_labels_20261013_scout.json
```

Artifacts: [fixed slices](../data/enumerations/fresh_h9_asymmetric_blocks_20261013_scout.json),
[three-label intersections](../data/enumerations/fresh_h9_asymmetric_intersections_20261013_scout.json),
[their selections](../data/enumerations/fresh_h9_asymmetric_intersections_20261013_selection_scout.json),
[four-label intersections](../data/enumerations/fresh_h9_asymmetric_four_labels_20261013_scout.json),
and [their selections](../data/enumerations/fresh_h9_asymmetric_four_labels_20261013_selection_scout.json).

This screen does not exclude caps with dotted seed pairings, more labels,
unsampled parameters, or an asymmetric two-by-two core metric with three
independent entries. In particular, choosing that third metric parameter
to enforce compatibility between cap families remains a distinct next
construction problem. No new-dimensional example follows from this scout.
