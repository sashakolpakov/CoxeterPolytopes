# A compact heptagon after releasing54 and replacing12

The exact 64-wall chamber improves the complete intrinsic S6 boundary from
17 to16 open ends. It has **197 finite vertices** and an actual compact
heptagonal two-face S6+8. A strict ambient escape survives, so this chamber
still has infinite volume.

This changes the [frozen release34 alternative](fresh_nonarithmetic_s6_release34_exact_20261017.md):
delete its wall54 and replace its wall12, retaining the other63 normals.
Wall labels retain their old meaning except12. The
[exact source](../scripts/fresh_nonarithmetic_release54_replace12_20261017.py)
and [certificate](fresh_nonarithmetic_release54_replace12_20261017.json)
give all native contacts, a common strict center, all64 actual-facet
witnesses, the complete finite-vertex census and the escape. The separate
[complete intrinsic boundary certificate](fresh_nonarithmetic_release54_replace12_boundary_20261017.md)
includes the costs of both deletions; it is not a full ambient exterior cone.

## Exact normal

Write S6={0,1,2,3,5,6}. Its white diagram is H4+A2. Old7 has the single
white product −1 at6, so its projection to S6-perp is

`e7=n7+(n5+2*n6)/3`, of square4/3.

Use the Lorentzian three-frame `(e7,n8,n52)`, with Gram

```
[ 4/3   -phi    e7.n52 ]
[ -phi    2        0   ]
[ e7.n52  0        2   ].
```

Let p be its uniquely determined vector with products `(0,0,-sqrt(3))`.
Let k be the old12 normal minus its projection onto this frame. Both k²
and2−p² are positive. The new norm2 normal is

`r=p+sqrt((2-p²)/k²)*k`.

It is right to S6,7,8 and makes angle pi/6 with52. Its remaining54 retained
native contacts are strictly dotted. The square2 and prescribed contacts
are identities from the orthogonal decomposition; interval overlaps are
only consistency checks. The chosen positive branch is strict against the
inherited timelike center. Projecting that center to each normal supplies
a finite point on every facet, strict against every other wall by acuteness.

Old54 is obtuse to r, so it really must be released in this recipe. Its
projection to S6-perp was future timelike, which explains the absence of an
intrinsic S6 deletion cost. Its ambient effect is checked separately below.

## Ambient finite vertices and preserved compact faces

The previous complete census had196 vertices. Old54 occurs at none of them.
Old12 occurs at exactly one:

`S6+[7,8,12]`.

The other195 old vertices retain their full positive-definite active sets.
The acute projection lemma makes the new inequality strict at every such
vertex. Removing old walls cannot reveal a further finite vertex made only
of retained old roots: every positive-definite full active set was already
an actual old vertex by the same lemma.

The new wall has exactly nine finite neighbors, S6+[7,8,52]. All nine
eight-element subsets are checked. Seven contain the dotted pair7,52 and
are impossible. The remaining two give actual finite vertices:

* S6+[7,8,new12];
* S6+[8,new12,52].

Thus the complete new finite census is195+2=197.

**No old compact face of dimension at least two is lost.** The deleted54
has no finite vertex. The sole lost old vertex has only one compact
incident edge, the edge supported on S6+[7,8]; every other incident edge
retains12 and cannot have a second finite endpoint, because no other old
vertex contains12. Consequently that vertex lies on no compact two-face
or higher compact face. The one affected compact edge receives the new
endpoint S6+[7,8,new12]. All old compact faces of dimension at least two
have their vertices preserved and strict against the new wall.

The inherited strict spacelike escape is future-oriented relative to the
common timelike center and remains strict against all64 walls. The open
ideal-boundary region it witnesses proves infinite volume. No construction
of a finite-volume H9 polytope is claimed.

Generation and fresh saved replay passed. The arithmetic uses the frozen
parent's outward rational interval engine and its isolated algebraic
parameter; this package is a new exact construction check, not an
independent implementation of that inherited arithmetic.
