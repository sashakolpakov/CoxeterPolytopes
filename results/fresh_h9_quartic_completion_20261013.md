# Full quartic anchors and two targeted H9 cap repairs

**Result:** 167 exact mirrors distinct from the preceding saved banks, and a
coherent progression from the existing 32-facet chamber to 35, 43, 44, 45 and
46 facets. Every tested chamber still has infinite volume, certified by a
strict future spacelike vector. No new finite-volume or compact H9 example
is claimed.

The last two additions are constructive repairs of specific open ends:
a new polar angle \(\pi/3\) closes saved target 3, and a seed angle \(\pi/4\)
closes the opposite gap at saved target 1. Both preserve every previous
facet. The final 46-facet chamber has a different exact escape, archived as
the next target. A subsequent cusp audit proves the stronger limitation
that **this fixed retained core has no finite-volume completion by adding
acute walls**. Continuing from its escape therefore requires replacing
obstructing retained walls. The underlying ten-wall seed remains open.

## Full quartic-anchor construction

The seed is the ten-wall, norm-two, signature \((9,1)\) seed obtained from
the Bugaenko H8 diagram by the earlier \(5\to3\) lift and deletion of
original wall 6. See the
[seed archive](../data/enumerations/fresh_h8_alternate_3_6_20261013.json).
Write \(K=\mathbf Q(\phi)\), \(\phi=(1+\sqrt5)/2\), and use dual coordinates
\(y_i=-(x,r_i)\). If \(G\) is the seed Gram matrix and \(I=G^{-1}\), then
the Lorentz norm is \(y^T I y\). A cap with seed pairings \(p_i\le0\) has
seed coefficients \(c=Ip\) and inequality \(c\cdot y\ge0\).

The retained 22-cap core is
`31,51,67,89,112,156,197,310,348,376,456,469,496,643,666,683,717,837,862,1028,1039,1067`.
Together with the ten seed walls this is the existing repaired 32-facet
configuration. Anchor 1039 has form \(A+\sqrt d\,B\), where
\(d=(-364+228\phi)/5\).

The new search solves the **full** equation
\(\langle r,A+\sqrt d\,B\rangle=-M\), together with a second anchor
equation. Three unspecified seed pairings leave a line over
\(E=K(\sqrt d)\); imposing \(r^2=2\) gives a quadratic over \(E\).
It does not impose the earlier additional condition \(\langle r,B\rangle=0\).
Targeted dotted distances \(M\) come from specific existing escape
directions. Subsequent rounds require compatibility with the actual
selected caps, so they preserve an increasingly large coherent chamber.
This is a proposal search, not an exhaustive enumeration of that family.

| Exact cap bank | Accepted entries | Required existing chamber | Tested enlarged chambers |
|---|---:|---|---|
| [targeted](../data/enumerations/fresh_h9_quartic_completion_20261013_targeted.json) | 80 | selected parts of the 32-facet core | six 35-facet chambers |
| [adaptive](../data/enumerations/fresh_h9_quartic_completion_20261013_adaptive.json) | 80 | all 35 facets of the first selected chamber | one 43-facet chamber |
| [second](../data/enumerations/fresh_h9_quartic_completion_20261013_second.json) | 5 | all 43 facets | one 44-facet chamber |

In the first 35-facet chamber, new caps 27, 74 and 76 are mutually dotted
and retain all 22 old caps. The adaptive 43-facet chamber retains this
entire chamber and adds eight further facets. The 44-facet chamber adds
second-bank cap 0.

The [exact mirror comparison](../data/enumerations/fresh_h9_quartic_completion_20261013_mirrors.json)
proves that all 165 entries are mutually distinct and distinct from all
2,865 preceding saved entries. The preceding count includes aliases;
it is not a count of 2,865 distinct mirrors. Equality is tested by all ten
normalized squared seed pairings. Nonpositive orientations and independence
of the ten seed normals make this a complete test of oriented mirror
equality. This is novelty relative to saved banks, not a worldwide
priority claim.

The accepted caps lie in quadratic towers of degree at most eight over
\(\mathbf Q\). Pair fields for these three banks have degree a power of two
at most sixteen. For a finite angle label \(m\ge3\), the degree
\(\varphi(m)/2\) of \(4\cos^2(\pi/m)\) must divide the pair-field degree;
therefore it lies in \(\{1,2,4,8,16\}\). The verifier checks every eligible
label, including the separately handled right angle:

`2,3,4,5,6,8,10,12,15,16,17,20,24,30,32,34,40,48,51,60,64,68,80,96,102,120`.

It verifies the polynomial identity
\(z^dP(2+z+z^{-1})=\Phi_m(z)\) and isolates the largest root by a Sturm
calculation, so a wrong conjugate is rejected. The finite scan to 2,048
is justified by \(\varphi(m)^2\ge m/2\) and \(\varphi(m)\le32\);
the power-of-two degree condition is an additional necessary filter.
The optional new-cap graph excludes new parallel pairs. Four parallel
relations are inherited from the retained core: old cap 31 with seed 8;
old cap 51 with seeds 8 and 9; and old cap 197 with seed 1. Each is a
norm-two pair with pairing \(-2\). These are valid Coxeter relations. In this
particular core, however, each forced cusp is also incident to the fixed
\(I_2(5)\) seed pair, which obstructs a finite-volume completion; see the
cusp audit below. The original infinite-volume certificates use spacelike
witnesses and do not rely on this later obstruction.

## A new polar-angle coupling: facet 45

The full quartic-anchor bank leaves saved target 3. A different coupling
releases the prescribed angle with anchor 1039 and keeps that relation
dotted. It sets only the two negative seed pairings

\[
w_4=\delta s,\qquad w_8=\delta(s+1),\qquad
\delta=\sqrt{2+\phi},\qquad p=-w.
\]

Let
\[
A=-20+\frac{25}{2}\phi,\quad B=-25+15\phi,\quad
C=-\frac{42}{5}+\frac{24}{5}\phi,\quad
D=B^2-4AC=-62+39\phi.
\]
Then
\[
s=\frac{-B+\sqrt D}{2A}
\]
is the unique positive root of \(As^2+Bs+C=0\):
\(A>0\), \(B<0\), \(C<0\), and exactly \(3<s<4\).
The resulting cap has norm two and angle \(\pi/3\) with polar facet 348.
It has eight right seed angles, two dotted seed relations, and dotted
relations to the other 33 retained nonseed facets. It cuts saved target 3
and has a strict common timelike interior with all 44 prior walls.

The [exact cap and 45-wall volume archive](../data/enumerations/fresh_h9_quartic_completion_20261013_polar.json)
proves that this mirror is distinct from the preceding 3,030 entries:
none has its nonzero seed-pairing support \(\{4,8\}\).
Cross-field pairings can have larger composita than the earlier
sixteen-degree bound. No finite-label bound is used for these new
cross-pairings: their negative signs and strictly dotted Gram determinants
are proved directly.

This chamber still admits a strict escape equal to saved targets 1 and 30.
Its underlying simpler two-coordinate ray is
\[
y_6=2-\phi,\qquad y_9=1,\qquad y_i=0\ (i\ne6,9),
\]
with norm \(1-\phi/2>0\) and center pairing \(-7/2\). It satisfies every
one of the 45 inequalities; old caps 67 and 310 have zero pairing.

## The opposite-gap repair: facet 46

The inverse seed Gram matrix has
\[
I_{66}=\frac{4+2\phi}{5},\quad
I_{99}=-\frac95+\frac{11}{10}\phi,\quad I_{69}=0.
\]
Set
\[
w_6=\sqrt2,\qquad
w_9=\sqrt{\frac{2-2I_{66}}{I_{99}}}=\sqrt{12+20\phi},
\qquad w_i=0\ (i\ne6,9).
\]
The cap \(c=-Iw\) has norm two, seed angle \(\pi/4\) at wall 6,
one dotted seed relation at wall 9, and eight right seed angles.
Every one of its 35 retained nonseed relations is strictly dotted.
It preserves the strict center and cuts both the displayed two-coordinate
ray and the strict saved targets 1 and 30.

The [exact opposite cap and 46-wall archive](../data/enumerations/fresh_h9_quartic_completion_20261013_opposite.json)
checks all 45 new pairings and mirror distinctness from the previous
3,031 entries. None has nonzero seed-pairing support \(\{6,9\}\).
This is the 167th new mirror of this round.

## Exact geometry and the remaining target

The final configuration consists of ten seed walls, all 22 old core caps,
twelve selected full-quartic-anchor caps, and the two explicit caps above.
The [44-wall geometry replay](../scripts/fresh_h9_quartic_completion_20261013_geometry_verify.py)
checks seed inertia \((9,1)\), the seed inverse, every seed and old/new
angle, all selected pairings, and the strict common timelike center.
The two extension verifiers check all additional relations.

These are true facets, not redundant inequalities. For any retained
spacelike normal \(r\), project the common strict timelike point \(x\)
orthogonally onto \(r^\perp\):
\(x_r=x-\langle x,r\rangle r/r^2\).
Its norm remains negative. Every other inequality remains strict because
all distinct outward-normal pairings are nonpositive. Thus the projected
point witnesses the relative interior of that facet.

Numerical halfspace vertices only propose escape vectors. Each archived
escape is then checked exactly: positive dual coordinates, positive Lorentz
norm, negative pairing with the timelike center, and strict satisfaction
of every cap inequality. Interpolating with the strict interior point
produces an open patch of feasible ideal directions, hence infinite
volume. No conclusion relies on numerical vertex completeness.

The final strict escape is rational:
\[
y=\frac1{10^9}
(22469,36355,36355,523224103,36355,36355,36355,36355,
408610843,68924455).
\]
Its exact norm is
\[
-\frac{10136812998738693821}{1250000000000000000}
+\frac{1579940765309177283}{312500000000000000}\phi>0,
\]
and its center pairing is
\[
-\frac{1846043559}{2000000000}
-\frac{477407883}{200000000}\phi<0.
\]
All 46 inequalities are strict. Exact cross multiplication proves that
this ray is not proportional to any of the earlier 36 saved target
entries. Its largest dual coordinates are 3, 8 and 9; that is a numerical
construction hint, not a claim that the strict vector has only three
nonzero coordinates.

The [next-target archive](../data/enumerations/fresh_h9_quartic_completion_20261013_next_targets.json)
preserves this witness as index 36, with its exact certificate and the
earlier 36 entries. It has 37 entries; inherited repetitions are retained.
No further construction pass was run on this new target in this round.
The cusp obstruction below means it is a target for a strategy that also
replaces retained walls, not a viable fixed-core wall-addition strategy.

## Finite volume requires changing the 46-wall core

The independently reconstructed
[four-wall cusp certificate and proof](fresh_h9_fixed_core_cusp_obstruction_20261013.md)
already applies to the retained 44-, 45- and 46-wall cores.

Each of the four inherited parallel pairs above consists of norm-two
normals \(a,b\) with \(\langle a,b\rangle=-2\). Its sum \(u=a+b\) is
nonzero and null, future-oriented by its negative pairing with the strict
timelike center. All retained and additional acute-compatible walls have
nonpositive pairing with \(u\), so the ideal point survives.

Exact seed pairings show that both seed normals 6 and 7 are orthogonal to
each such \(u\). They have Gram matrix
\[
\begin{pmatrix}2&-\phi\\-\phi&2\end{pmatrix},
\]
so the induced Euclidean cusp link retains their angle \(\pi/5\).
For any proposed full-dimensional completion, its own strict timelike
interior point induces a strict interior point of the link. The archived
center certifies the original null ray's orientation and nonvanishing;
newly added walls need not retain that particular center.
Orthogonal projection of the completion's link point onto either incident wall
preserves all other inequalities strictly; simultaneous projection onto
these two walls does the same because their Gram inverse is nonnegative.
Thus these are true link facets meeting at a true ridge.

A finite-volume completion would require a bounded Euclidean
eight-dimensional cusp link. In upper-halfspace coordinates with this
cusp at infinity, incident walls are vertical and the finitely many
nonincident hemispheres disappear above their maximum radius. An
unbounded vertical-wall link would therefore give infinite volume.
But a bounded Euclidean Coxeter chamber yields a crystallographic
reflection group. Its translations form a lattice preserved by the
linear parts of its reflections, by Bieberbach's theorem; see
[Gelander, Lectures on Lattices, Theorem 2.9](https://www.math.utah.edu/pcmi12/lecture_notes/gelander.pdf).
The product of reflections in the two
\(\pi/5\) walls has linear trace
\[
6+2\cos(2\pi/5)=5+\phi,
\]
which is irrational, whereas a lattice automorphism has integer trace.
This contradiction rules out every finite-volume acute wall-addition
completion retaining this core. It does not rule out replacing walls or
completing the ten-wall seed by a different core.

## Deleting parallel pairs still leaves compactness obstructions

The four inherited parallel pairs already prove that the fixed 46-wall
core cannot be made compact by adding acute walls. If a set of retained
normals has positive coefficients \(\lambda_i\) and a singular
positive-semidefinite Gram matrix with kernel \(\lambda\), then
\(v=\sum_i\lambda_i r_i\) is nonzero and null. Its pairing with the strict
timelike center is negative, and its pairing with every retained or new
acute-compatible normal is nonpositive. The feasible ideal direction
therefore survives every such wall addition.

The [compactness-scope archive](../data/enumerations/fresh_h9_quartic_completion_20261013_compact_scope.json)
removes old caps 31, 51 and 197. This gives a specified 43-facet
parallel-free Coxeter chamber with the same strict spacelike escape.
It is a separate deletion checkpoint, not the earlier adaptive 43-facet
chamber. Five positive affine circuits remain:

| Type | Retained walls |
|---|---|
| \(\widetilde A_3\) | seed 0, seed 1, seed 2, adaptive cap 19 |
| \(\widetilde A_3\) | seed 2, seed 3, seed 4, old cap 1028 |
| \(\widetilde A_3\) | seed 5, seed 6, old caps 310 and 837 |
| \(\widetilde A_2\) | seed 3, seed 4, old cap 496 |
| \(\widetilde A_2\) | seed 4, seed 5, old cap 496 |

Every displayed Gram matrix is checked exactly, with positive kernel
consisting entirely of ones. Hence this 43-wall core also cannot be
made compact by acute wall addition. With the ten seed walls fixed,
the smallest sets hitting these five displayed circuits consist of
adaptive cap 19, old caps 1028 and 496, and one of old caps 310 or 837.
This is a minimum only for the displayed circuits; removing such a set
is not a certified compact construction and may leave other affine
circuits or open ends.

The ten-wall seed itself has no affine principal subdiagram. Its
8–9 principal determinant is negative, so a positive-semidefinite subset
must omit one of those walls. The two nine-wall principal matrices
obtained by omitting 8 or 9 are positive definite, checked by exact
elimination (types \(A_7+I_2(5)\) and \(E_7+I_2(5)\)). Every remaining
principal subset is therefore positive definite. The five displayed
circuits are not orthogonal to both seed normals 6 and 7, so the preceding
\(I_2(5)\) cusp argument does not automatically transfer to this 43-wall
deletion checkpoint. No finite-volume nonexistence claim is made for
that checkpoint. The compact target requires further changes to retained
caps or seed walls.

## Replay and independent review

All commands run from the repository root with `PYTHONPATH=.`. The first
three accepted-root archives have producer-free Fraction replays:

~~~sh
python scripts/fresh_h9_quartic_completion_20261013_verify.py --verify
python scripts/fresh_h9_quartic_completion_20261013_verify.py --source data/enumerations/fresh_h9_quartic_completion_20261013_adaptive.json --verify
python scripts/fresh_h9_quartic_completion_20261013_verify.py --source data/enumerations/fresh_h9_quartic_completion_20261013_second.json --verify
python scripts/fresh_h9_quartic_completion_20261013_mirrors.py --verify
python scripts/fresh_h9_quartic_completion_20261013_geometry_verify.py --verify
python scripts/fresh_h9_quartic_completion_20261013_geometry_verify.py --source data/enumerations/fresh_h9_quartic_completion_20261013_adaptive_geometry.json --verify
python scripts/fresh_h9_quartic_completion_20261013_geometry_verify.py --source data/enumerations/fresh_h9_quartic_completion_20261013_second_geometry.json --verify
python scripts/fresh_h9_quartic_completion_20261013_polar.py --verify
python scripts/fresh_h9_quartic_completion_20261013_opposite.py --verify
python scripts/fresh_h9_quartic_completion_20261013_compact_scope.py --verify
python scripts/fresh_h9_fixed_core_cusp_obstruction_20261013.py --verify
~~~

The polar, opposite-cap and compactness-scope `--verify` paths import no
numerical libraries or search producers. A second agent independently reconstructed both cap formulas,
all new pairings and both full-volume witnesses using a separate
compositum implementation, and rejected corrupted cap/witness data:
[polar independent audit](fresh_h9_quartic_completion_20261013_polar_independent.md),
[opposite independent audit](fresh_h9_quartic_completion_20261013_opposite_independent.md).
A third reviewer checked the tower arithmetic and independently derived
the formulas:
[formula review](fresh_h9_changed_form_20261013_polar_review.json).

The evidence concerns the specified exact configurations and saved
mirrors. The fixed 46-wall core cannot be completed to finite volume by
acute wall addition. A finite-volume completion of the ten-wall seed
after suitable changes to the retained cap core remains unexcluded.
