# An exact H20 cap with a short finite chain arc

The frozen133-wall chamber admits a coherent134th wall. It has angles
π/4,π/3,π/4 with chain walls128,129,130, is right to17 support walls,
and is dotted to all113 remaining walls. It cuts the strict escape previously
certified against the847-recipe nonright-support bank. It adds exactly two
finite vertices on its facet and three finite edges connecting them to two
old finite vertices. The chamber remains infinite volume.

Primary generation and fresh saved replay pass. This is a new construction
certificate; an independent arithmetic review is pending at publication of
this report. The numerical scout was used only to select the literal recipe.

Artifacts: [exact checker](../scripts/fresh_h20_dotted_arc_cap_20261016.py),
[exact result](fresh_h20_dotted_arc_cap_20261016.json). Dependencies include
the [133-wall construction](fresh_h20_s61_multicap_20261014.md), its
[independent review](fresh_h20_s61_multicap_20261014_exact_review.md), and the
[847-recipe escape certificate](fresh_h20_nonright_s61_20261014.md).

## Exact definition

Use the frozen orthonormal Lorentz frame `(x,e,f)` perpendicular to the
rank18 support S61. Let `a` be the existing isolated positive algebraic
parameter, and write the central chain normal as `n129=(a,b129,c129)`.
Set

```
h = sqrt(2)-1
y = ((h*a^2+2*sqrt(2))/(2*a), h*b129/2, h*c129/2)
c = sqrt(5+4/a^2)
r = y-c*omega31,                 omega31=n31/2.
```

Root31 is an isolated A1 component of S61, so `omega31²=1/2` and this
dual vector is perpendicular to the Lorentz frame and all other support
roots. The saved rotation has cosine `C=a²/(a²+2)`. Its adjacent spatial
dot products immediately give

```
y.n129 = -1,             y.n128 = y.n130 = -sqrt(2),
y² = -1/2-2/a²,          r² = y²+c²/2 = 2.
```

These are algebraic identities, not acceptance from intervals containing
zero. Every unprescribed pair is certified strictly less than−2 using
448-bit outward rational intervals. The largest dotted product is less
than−2.10336. All26 chain walls are retained: three finite contacts and23
dotted contacts. The nonabelian-chain obstruction is therefore respected.

## Coherence and actual facets

The old strict center `H=x-e/2+epsilon*H107/sqrt(-H107²)` remains strict.
Its new pairing is approximately−2.81319865. All134 center pairings and
the center's negative norm are certified.

For any normalized wall `n_i`, the vector
`H_i=H-(H.n_i)n_i/2` is future timelike and orthogonal to `n_i`.
For every other wall `n_j`, acuteness gives
`H_i.n_j=H.n_j-(H.n_i)(n_i.n_j)/2 < 0`.
Thus all134 walls define actual facets; none is redundant. The same acute
extension argument preserves the168 previously certified finite vertices.
This is a lower bound on the old ambient census, not a claim of completeness.

## The two new finite vertices and three new edges

The complete old finite-neighbor list of the new root is
`(S61 minus31) union {128,129,130}`. A finite vertex on its facet needs19
old incident roots. Every19-subset except omission128 or omission130
contains the dotted pair128/130 and is nonspherical. For each of the two
remaining supports, the last three-root block has bonds√2 and1 to the
new root and zero between the old roots; its Schur remainder is1/2.

Let `z` be the old S61 polygon vertex at the adjacent pair128/129 or129/130,
normalized to time coordinate1. The new vertex is explicitly

```
v = z + (2*(z.y)/c)*omega31.
```

Both points have squared norm approximately−0.000706898253 and are strictly
negative on every nonincident retained wall. The defining orthogonality
relations are exact. This proves exactly two finite vertices on the new
facet and at least170 finite vertices overall.

The new vertices share17 support roots, wall129 and the new wall, giving
their finite edge. Each is joined to the corresponding old polygon vertex
by the support consisting of17 retained support roots and its adjacent
chain pair. All endpoints are feasible future timelike; the open segments
are strict on the remaining walls. These are three distinct new finite
edges. Their four endpoints form a compact quadrilateral section at
`(S61 minus31) union {129}`. No complete ambient edge census is asserted.

## Escape status

The old847-bank escape is cut with pairing greater than0.49838. The still
older133-wall escape survives with pairing less than−32.0327. All three
previously certified proper cusps also remain strict.

A simpler additional strict escape is included to give a self-contained
volume witness. In the same orthogonal frame and support-dual basis, take

```
Enew = x-e/100-omega5-(1/100000)*sum(omega_j : j in S61, j != 5).
```

Its norm is the exact positive rational saved in the JSON (approximately
1.72761819); all134 wall pairings are strictly negative. Its pairing with
the new cap is less than−3.54196, and its center pairing is negative. Thus
the134-wall chamber still has infinite hyperbolic volume.

The result is a concrete escape-cutting extension with new finite geometry.
It supplies no finite-volume polytope, new dimension, or all-parameter
classification of short-arc caps. The earlier450-family scout and the
bounded triple-label search are exploratory discovery data and are not
dependencies of this exact proof.
