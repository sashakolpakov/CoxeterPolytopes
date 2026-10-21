# Exact closure of a three-cap patch after changing two old walls

This alternative deletes old walls 34 and 52 from the frozen 63-wall
chamber, retains the other 61, and adds a replacement 52 together with
new walls A=63, B=64 and C=65. There are **65 actual facets**; the labels
are not consecutive. Three entire two-faces are compact: two quadrilaterals
and a pentagon. A strict remaining escape proves that the ambient H9
chamber still has infinite volume.

The [exact producer](../scripts/fresh_nonarithmetic_s6_release34_exact_20261017.py)
and [certificate](fresh_nonarithmetic_s6_release34_exact_20261017.json)
use literal algebraic definitions and outward rational intervals. The
[numerical scout](fresh_nonarithmetic_s6_release34_patch_scout_20261017.md)
motivated the construction but supplies no exact acceptance test.

## The released parameter and its exact definition

Put `S6=[0,1,2,3,5,6]`, of type H4+A2, and write
`theta_m=2cos(pi/m)`. All four new normals are right to S6. Their other
prescribed finite products are:

| Normal | Right contacts outside S6 | Nonright finite contacts |
|---|---|---|
| Replacement 52 | 8 | 30: pi/4; C: pi/4 |
| A=63 | 30,44,C | B: pi/27 |
| B=64 | 30,replacement52,C | A: pi/27 |
| C=65 | 8,A,B | 44: pi/15; replacement52: pi/4 |

The contact of replacement 52 with B is also right. Every unspecified
contact with a retained or new normal is strictly dotted.

Let `b=<old8,old30>` and let x be the free dotted product
`<old44,replacement52>`. For the ordered frame `[8,44,replacement52]`, set

```
G3 = [[2,0,0], [0,2,x], [0,x,2]],
c = [0,-theta15,-theta4],
d = [b,0,-theta4],
rho = 2-c^T G3^-1 c,   w2 = 2-d^T G3^-1 d,
z = c^T G3^-1 d - sqrt(rho*w2).
```

The Gram on `[30,44,replacement52,C]` is then

```
G4 = [[2,0,-theta4,z],
      [0,2,x,-theta15],
      [-theta4,x,2,-theta4],
      [z,-theta15,-theta4,2]].
```

With `J=G4^-1`, define x as the unique root of

`2 J[1,2]/sqrt(J[1,1] J[2,2]) + theta27 = 0`

in the rational interval `[-4.098064174,-4.098064173]`. The checker proves
opposite signs at its endpoints and a derivative throughout the interval
between approximately 5.2781998958 and 5.2782003943. Every inverse pivot
and positive radicand is certified over the whole interval. Forty-four
rational bisections isolate the root to width less than `1e-22`, near
`-4.0980641738104895`.

The equation is an algebraic function of x and the fixed algebraic old
data. The sign and derivative certificates establish existence and
uniqueness of an exact algebraic root. The root is not defined by rounding
the numerical scout, and an interval merely containing a Coxeter value is
not used as an equality proof.

## Reconstructing the four normals

The three old normals `[8,30,44]` have nonsingular Gram
`[[2,b,0],[b,2,0],[0,0,2]]`. Their prescribed products with replacement 52
are `[0,-theta4,x]`. Let p(x) be the corresponding projection and let k
be the component of old 52 perpendicular to this frame. Both k and p(x)
are perpendicular to S6. The checker proves `k^2>0` and defines

`replacement52 = p(x) + sqrt((2-p(x)^2)/k^2) k`.

Using the removed old 52 as an algebraic coordinate vector does not retain
its old halfspace. The positive branch agrees with that old vector at the
old value of x.

Next project the prescribed contact row c onto `[8,44,replacement52]`.
Subtract the positive normalized component of old 30 perpendicular to this
frame to obtain C; its product with old 30 is exactly z above. Finally,
A and B are the negative normalized dual vectors to replacement 52 and
44 in G4. Orthogonal projection and inverse-Gram identities prove their
norms and all prescribed finite contacts. The defining equation for x
proves the remaining A--B angle is exactly pi/27.

The certificate checks **250 new contacts**: 214 strictly dotted, 32 right,
two pi/4, one pi/15 and one pi/27. Old/old contacts are inherited after the
two explicit deletions. The old future timelike center remains strict.
For every normal n, its orthogonal projection of that center onto `n^perp`
is timelike and strict at all other walls by acuteness. Thus all 65 listed
walls are actual facets.

## Complete polygons and finite vertices

Adjoin S6 to each of the following eight side triples:

```
[8,44,C], [8,replacement52,C],
[30,44,A], [30,replacement52,B],
[44,A,C], [replacement52,B,C], [30,A,B], [A,B,C].
```

Each resulting nine-root Gram is positive definite, and its future
timelike vertex is feasible with every nonincident wall strict. Norms
and wall products are checked by exact projection identities, avoiding
cancellation in expanded physical coordinates.

The complete side cycles of the three new two-faces are:

| Face support | Cyclic sides |
|---|---|
| S6+A | 30,44,C,B |
| S6+B | 30,replacement52,C,A |
| S6+C | 8,44,A,B,replacement52 |

For each support, exactly these side normals have positive projected
square. Every other projected wall normal is timelike, and its halfspace
contains the whole relevant future hyperbolic plane: acute projection of
the strict center fixes the orientation. The listed cyclic corners are
finite and strict at every other side. Consequently each side segment is
bounded at both ends, and their convex polygon is the entire face. This
certifies whole quadrilaterals and a pentagon, not isolated corner diagrams.

The complete ambient finite-vertex count is **196**. Exactly five of the
old 193 finite corners use deleted 34 or old 52. The other 188 are checked
strict against all four new normals. Exhausting every eight-neighbor
subset at each new wall yields precisely the eight listed new corners;
every other candidate contains a certified dotted pair. A finite corner
using only retained old walls already belonged to the old acute chamber,
by projection of its strict center onto the positive-definite nine-root
support. Thus no further old-only finite corners appear after deletion.

## Remaining scope and verification

The old 63-wall strict escape remains spacelike, future relative to the
center, and strict against all 65 inequalities. Its squared norm is about
0.01400393656. The segment from the center to this escape crosses the
ideal boundary strictly inside every halfspace, proving infinite volume.
There is no claim of a complete ambient exterior cone or persistence of
old faces involving either deleted wall.

The separate [numerical boundary ledger](fresh_nonarithmetic_s6_release34_boundary_20261017.md)
counts 17 open S6 ends, down from 19, while accounting for the new opening
where 34 was removed. An exact boundary census is a separate successor.

Validation status: **FINAL/FROZEN**. Generation and a fresh saved replay passed,
including the complete finite-vertex census. A separate agent's narrow mathematical/source review
passed for the algebraic-root equality and whole-polygon arguments. That
review is not an independent arithmetic reconstruction or external refereeing.
