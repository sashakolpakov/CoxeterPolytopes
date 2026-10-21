# Changed weighted-box seed and two covering obstructions

The new seed uses nineteen affine `A1` factors with physical inradii
`1` (ten factors), `sqrt(2)` (one), and `phi` (eight), where
`phi=(1+sqrt(5))/2`. Its 38 vertical walls and the radius-two central
hemisphere form a coherent 39-wall H20 Coxeter core, with central labels
3, 4, and 5. This changes the cusp factors as well as their relative scales.
It does not produce a finite-volume example.

The exact [certificate](fresh_h20_seed_deformation_20261014.json) and
[standalone producer/replay](../scripts/fresh_h20_seed_deformation_20261014.py)
establish two **separate** global obstructions. They do not exclude mixing
the two families, other radii, seed labels above 6 in the second family,
or changing the central cap. Generation and saved `--verify` both passed.

## Coordinates and the scale gate

Write physical coordinates as `b_i x_i`, so the old Euclidean box is
`[-1,1]^19` and its spatial Gram matrix is diagonal with entries
`1` ten times, `2` once, and `phi^2=1+phi` eight times. Every box vertex
has square `N=20+8phi`. In the Lorentz model use

\[
 (x,d,u)^2=x^2-2du.
\]

The vertical normals are `(±e_i,b_i^2,0)`, and a sphere of center `c`
and radius `R` has normal `(c,(c^2-R^2)/2,1)`, square `R^2`.
The central cap is `(0,-2,1)` and the common strict timelike point is
`(0,100,1)`. These formulas give the stated signature and Coxeter angles
directly.

There is a useful restriction on deforming the earlier product seeds.
If a central sphere of radius `R` has its center strictly inside every
seed inequality, every positive Coxeter seed distance is at least `R/2`:
the smallest nonright finite cosine is `cos(pi/3)=1/2`, and a dotted
distance is larger. Taking the positive affine weighted average of the
distances in each factor gives `b_f >= R/2`. For an `A`-factor product,
the vertices have common norm and barycenter zero. Their average squared
distance from an off-center point `c` is therefore `N+|c|^2`. Thus moving
an interior central center cannot improve the worst vertex-distance
ratio below the centered unit-scale value for the same factors.
This argument does not cover boundary-centered caps with zero seed gaps.

## Radius-four vertex family

Put a radius-four hemisphere at any of the `2^19` box vertices. Each is
right to nineteen seed walls, has labels 3/4/5 against their opposite
walls, and has label 5 against the central cap because

\[
 \frac{N-4-16}{2\cdot2\cdot4}=\frac{\phi}{2}.
\]

For two signed vertices let `h_1,h_2,h_5` count the differing unit,
sqrt-two, and golden coordinates. Their distance squared is `4H`, where
`H=h_1+2h_2+(1+phi)h_5`. The negative normalized cap pairing is
`(H-8)/8`. It lies in `Q(phi)`, whose possible finite Coxeter cosines are
only `0,1/2,phi/2`: the degree formula for `2cos(pi/m)` first leaves
`m=2,3,4,5,6`, and labels 4 and 6 require different quadratic fields.
Parallel/dotted relations require `H>=16`.

Within one golden sign pattern, compatibility forces weighted Hamming
distance 8 or 12. If there are `M` such caps, the weighted Plotkin sum
bound gives

\[
 8\binom M2\le 12M^2/4,
 \qquad M\le4.
\]

Golden sign patterns differing in one bit are incompatible for every
choice of the other coordinates. A matching in the golden 8-cube thus
bounds the number of represented golden patterns by 128, and any
coherent selection has at most 512 caps.

Each closed cap ball contains exactly 530 old box vertices, determined
by `H<=4`. Consequently at most `512*530=271360` of the 524288 box
vertices can be covered. The central ball covers none. A missed vertex
is strictly outside every selected closed ball; a sufficiently small
inward displacement in the box remains uncovered. This produces an
open ideal-boundary patch and hence infinite volume. The result covers
every coherent subset of this complete signed radius-four family.

## Radius-two nonvertex family

Now permit arbitrary radius-two centers compatible with all seed walls,
using seed finite labels 2 through 6, plus parallel/dotted relations.
The two opposite-wall gaps sum to `2b_i`. Exact one-coordinate solving
gives the following physical coordinates:

| Inradius | Allowed center coordinates |
|---|---|
| `1` | `0, ±1` |
| `sqrt(2)` | `0, ±sqrt(2)` |
| `phi` | `0, ±(phi-1), ±phi` |

No bound on the central label is imposed. A radius-two center's norm is
in `Q(phi)`, so its finite central cosines and its equal-radius cap-cap
cosines are again only `0,1/2,phi/2`.

Consider a half-vertex probe: in normalized coordinates each coordinate
is `±1/2`. Start with all golden center coordinates equal to the small
matching value `±(phi-1)`. The minimal squared distance is `13-6phi<4`.
Replacing a golden coordinate by zero or its full matching value adds
`phi-1`; two replacements already make the distance exceed 4. Any
opposite nonzero sign adds at least 2, so cannot cover the probe.

The all-small and one-zero golden cases fail central compatibility:
their maximum center norms are respectively `28-8phi` and `26-7phi`,
both below 16, and their candidate finite central cosines have a negative
golden coefficient. In the remaining one-full-golden case central
compatibility requires center norm at least 16. This forces the heavy
coordinate to be full and at most one of the ten unit coordinates to be
zero. Thus every core-compatible cap covering the probe is one of exactly
88 recipes `(a,i)`: one full golden coordinate `a`, seven small golden
coordinates, and either no zero unit coordinate or a single zero `i`.
Every recipe has strictly negative probe power `8-5phi`; there are no
boundary-only covering exceptions.

Two probes differing only in golden sign `j` have compatible cutters
exactly when one cutter, but not both, has its full golden coordinate
at `j`, and their two zero-unit coordinates are distinct and present.
The 7744 exact cross-pair checks give 1260 right pairs and 6484 forbidden
pairs. A cutter cannot cover both probes.

Take the eight probes obtained by varying just three golden signs.
Their adjacency graph is a 3-cube with twelve edges. Choose a covering
cap for each probe. Every adjacent pair requires at least one endpoint's
full-golden coordinate to designate that edge direction. Each cap can
designate at most one incident cube edge. Eight caps cannot supply twelve
edges. Therefore a coherent family cannot cover all eight probes.
The missed probe lies strictly inside the old box and strictly outside
all selected balls, again proving infinite volume.

The certificate includes a literal local cutter: one full golden
coordinate, one zero unit coordinate, and all other coordinates as above.
It is right to 11 seed walls, has nine label-3 seed contacts, nine
parallel and nine dotted seed relations, and is dotted to the central
cap. It cuts its half-vertex probe. The obstruction concerns joint
covering, not the absence of individual useful caps.

## Next scope

These tests show why merely changing balanced simplex factors to a
smaller weighted box does not close either of the two natural fixed-radius
families. A mixed-radius construction is still a distinct possible route:
the radius-four caps can cover probes not assigned a radius-two cutter,
so the eight-probe counting proof must not be applied to their union.
No lattice-root enumeration, arithmetic existence assertion, or global
nonexistence theorem is used.
