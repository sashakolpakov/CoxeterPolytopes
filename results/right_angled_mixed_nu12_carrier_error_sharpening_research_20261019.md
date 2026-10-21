# Mixed RA carrier errors: coefficients one and two without finite-edge terms

19 October2026. **New written successor, awaiting independent review.**
No computation, model or graph census, optimizer, Lean or build ran.
Existing source/certificate bytes are unchanged. This is a necessary
incidence theorem; the unrestricted mixed dimension ceiling remains12.

For an ambient finite-volume right-angled n-polytope let

    E12=sum_(three-faces X with nu(X)>=12)(nu(X)-11),

where nu(X)=F(X)+c(X). At a finite vertex let M count unbacked induced
two-edge paths in its bad-triangle graph. Let A3(X) count finite vertices
of X all of whose three incident edges have ideal other endpoint, and
omit the canonical P3 and its reflected double D3 from the A3 sum.
Then the strengthened universal estimates are

    sum_q M(q)<=E12,                            (1)
    sum_(X not P3,D3) A3(X)<=2E12.              (2)

No finite--finite edge term occurs. No global restriction on edges or
three-face types is imposed.

## 1. Exact inherited local inputs

The immediately preceding written theorem is
`right_angled_mixed_nu11_all_ideal_corner_exclusion_research_20261019.md`,
SHA256
`fc047b78e9a8fc9eaf283bd5d689b69d72a8bbf0976a4765600ffbe8625955b0`.
Its three-unit polygon-excess proof gives: a non-P3 right-angled
three-polytope with an all-ideal-edge finite corner satisfies
F+2c>=15; hence such a corner at nu<=11 occurs only in D3. Consequently
an unbacked induced-path corner, or an all-ideal-edge corner outside
P3,D3, can occur only at nu>=12. This is a LOCAL theorem and has no
global no-finite-edge premise.

Other inherited inputs are finite link degree3, cusp degree4, simplicity
of the edge graph, planar graph topology, and the actual attachment
lemma: two bad triangles sharing their ideal--ideal edge force P3.
A bad triangle means one finite and two ideal vertices. These are
written geometric inputs, not consequences of a numerical model check.

For every three-polytope Euler gives

    p+4c=2nu-4,
    4c=E_FI+2E_II,
    3p=E_FI+2E_FF.                             (3)

All edges here are actual distinct edges of its simple convex graph.

## 2. At most one induced-path corner at nu12

Let Z(X) count its induced-path finite corners. Every such corner has
three distinct ideal neighbors and lies in two bad triangular facets.
Different bad triangles use different ideal--ideal edges unless X=P3.
Therefore

    E_FI>=3Z, E_II>=2Z, 4c>=7Z, p>=Z,
    4Z<=nu-2.                                  (4)

Suppose nu12 and Z>=2. Equations(3)--(4) force c4,p4,F8: c>=4
from4c>=14, while p>=2 and p+4c20 force c<=4.
Choose two path vertices q,r. Both have all three incident edges ideal,
and together supply six finite/ideal edges. Their four bad triangles
use four distinct ideal/ideal edges, so(3) gives

    E_FI<=16-2*4=8.

But the remaining two finite vertices can have no finite-edge neighbor
among q,r. Each can have at most the other remaining vertex as a finite
neighbor. Each therefore has at least two ideal neighbors. Thus

    E_FI>=3+3+2+2=10,

a contradiction. Hence **Z<=1 at nu12**. This excludes a corner pattern,
not the entire F8,c4,p4 class or all of its possible finite-edge graphs.

At nu13, the integer form of(4) gives Z<=2=nu-11. At nu>=14,
(nu-2)/4<=nu-11 gives the same inequality directly. Thus

    Z(X)<=nu(X)-11                 for all nu>=12.       (5)

An unbacked path is among these Z corners; its nu<=11 contributions
are absent by Section1. Summing(5) proves(1). The coefficient2 in the
earlier error bound is no longer needed, even on nu12 faces.

## 3. All-ideal-edge corners at nu12 and nu13

Write A=A3(X). Every such finite vertex has three distinct cusp
neighbors. In particular c>=3 when A>0, and

    3A<=4c, p>=A,
    A<=(nu-2)/2.                               (6)

We now check the two small exceptional values needed for coefficient2.

At **nu12**, suppose A>=3. Then c>=3 and p>=3 with p+4c20 force
c3 or4. If c3, any three of these finite vertices are adjacent to
all three cusps, giving a K3,3 subgraph of the planar edge graph.
This is impossible.

If c4, then p4. Any fourth finite vertex cannot have a finite edge
to one of the first three, and has no other finite candidate. It too
has three ideal neighbors. Thus all four finite vertices have only
ideal-edge neighbors, E_FI12 and E_II2.

The facet corner total is3p+4c28 across F8 polygonal facets. Since
every nontriangle has at least four vertices, there must be at least
4F-28=4 triangular facets. A triangle with two finite vertices needs
a finite--finite edge, absent here. An all-ideal triangle needs three
distinct ideal--ideal edges, but there are only two. Thus every triangle
would be bad. The attachment lemma allows at most one bad triangle
per ideal--ideal edge, hence at most two, a contradiction. Consequently

    A<=2=2(nu-11)                         at nu12.       (7)

At **nu13**, suppose A>=5. Now p+4c22 and p>=5 again force c3 or4.
The c3 case is ruled out by the same K3,3 subgraph. If c4 then p6.
The sixth finite vertex has no finite-edge neighbor among the first
five, and is therefore also all-ideal-edge. This would give E_FI18,
exceeding the cusp capacity4c16. Hence

    A<=4=2(nu-11)                         at nu13.       (8)

For every nu>=14, (6) directly gives A<=2(nu-11). There is no relevant
nu<=11 contribution by Section1. Summing(7)--(8) and the general tail
proves(2).

## 4. Scope of the gain

The universal quad displacement still uses actual symmetric quad
moves weighted by half their invariant common-bad-neighbor multiplicity.
Several carriers sharing one missing pair do not become different
opposite vertices. Equations(1)--(2) alter the error ALLOWANCE, not
the definition or reversal multiplicity of these actual moves.

The comparison is:

|Error|Earlier reviewed allowance|New allowance|
|---|---|---|
|Unbacked paths|2E12+binom(n-1,2)e_FF/5|E12|
|Noncanonical all-ideal-edge triple corners|5E12+(2/3)binom(n-1,2)e_FF|2E12|

The previous nu11 error class was a conservative possibility, not an
asserted geometric realization. It is now excluded at these particular
corners by a polygon-excess argument. Other nu11 finite-edge carriers
are not excluded.

These sharper inequalities must be composed with a quantitative causal
or reciprocal argument and a complete ordinary count before excluding
dimension12. A single height-maximizing corner does not control the
unweighted global average; no such replacement is made here. Graphs
with a spacelike hypothetical full field also remain an independent
obstruction. The new ideal5/6 size bounds and ideal-face fractions
are unchanged. No unconditional dimension improvement follows from
this report alone.
