# A dyadic D4 extension with a complete first F4 edge

20 October 2026. Own paper construction and finite-interval proof.
No worker, height search, numerical calculation, Lean or TeX job.
This supplies a genuine actual edge, not a compact H9 example.

Put `phi=(1+sqrt(5))/2`, `O=Z[phi]`, and

    R_n=<-phi> + I_n + D4(O),                  n>=5,
    D4(O)={z in O^4 : z1+z2+z3+z4 in2O},

with ordinary dot product on the displayed D4 coordinates. The target
R_5 has hyperbolic dimension9. Its determinant is `-4phi`, so this
does NOT claim a new K-valued determinant squareclass:4 is square.
It is a different integral lattice of index ideal(2) in the odd
unimodular diagonal lattice. The quotient has cardinality4 as a
Z-lattice. No integral equivalence, noncommensurability, or worldwide
novelty claim is inferred from this index.

## 1. Complete norm pool and zero-height arrangement

The elementary golden arithmetic used here is O PID, units `+-phi^j`,
totally-positive units squares, and inert2. The dual exponent of R_n
is2: D4 has exponent2 and the other summands are unimodular. Primitive
reflection integrality therefore gives the complete necessary norm
pool, after unit rescaling,

    q in {1,2,4}.                                  (1)

For q4, the time and I_n coordinates lie in2O, and the D4 coordinate
lies in `2D4^* tensor O`. The ordinary lattice `2D4^*` has minimum4.
At time zero, nonzero I_n coordinates in2O already contribute trace
at least8, the entire trace of q4, leaving only the imprimitive
coordinate vectors2e_i. In the D4 block write z=a+phi b. Rational
norm q gives

    a^2+b^2=q,                 2(a,b)+b^2=0.        (2)

At q4, a,b lie in2D4^*, so (2) forces b=0 and a^2=4. These are the
ordinary24 vectors `+-2e_i` and `(+-1,+-1,+-1,+-1)`. Norm2 D4 roots
are the ordinary24 vectors `+-e_i+-e_j`. Norm1 roots occur only in
I_n; its norm2 roots are the ordinary B_n directions. Equations(2)
also exclude every mixed zero root. Thus the FULL zero chamber at
`p=(1;0)` has finite diagram `B_n+F4`.

Choose its F4 simples

    a1=e1-e2, a2=e2-e3, a3=2e3, a4=-e1-e2-e3+e4,

of norms2,2,4,4 and unit diagram(3,4,3). All four belong primitively
to D4. In particular the norm4 spinor wall reflects a coordinate of
the odd unimodular overlattice to a half-integral vector. Its mirror
arrangement cannot simply be identified with that overlattice's
norm1/2 arrangement. The known exclusions of specified unimodular
hosts do not establish nonreflectivity of R_n by themselves.

The F4 face has the FULL `<-phi>+I_n` arrangement. No connected
spherical extension of F4 exists, so its actual neighbor walls are
orthogonal to the whole F4. Norm1/2 intrinsic roots extend, while a
norm4 root orthogonal to D4 has every coordinate divisible by2 and
is imprimitive. Thus the design uses the known compact dimension5
diagonal face when n=5; it does not require a metric Bugaenko8 facet.
Compactness of this smaller diagonal face is a published input, not
a claim that compactness of the extension follows.

## 2. Explicit finite upper

Release a4, retaining B_n and the B3 on a1,a2,a3. The edge ray is

    P(t)=(1;0^n;0,0,0,-t).

The primitive integral reflector

    r=(2;0^n;0,0,0,-2phi),            r^2=4,

is orthogonal to every retained wall and has `r.a4=-2phi`.
It meets the ray at t=1, strictly before light t=sqrt(phi), with
`P(1)^2=1-phi<0`. Reflection integrality follows directly: the
unimodular coordinates are divisible by2 and all pairings of its
D4 part with D4 lie in2O. It is primitive because division by2
would give a D4 coordinate with odd sum `-phi`.

## 3. Complete finite first-wall profiles

Local finiteness gives an actual first wall no later than r. Orient
it to have positive time x. The retained whole components are B_n
and B3. A finite extension can connect to at most one of them: a
connected finite diagram cannot contain both of their4-bonds.
The finite classification allows only B_(n+1) from B_n, and B4 or
F4 from B3. Multiple contacts within a component create a cycle.
The resulting full inventory is:

|type|q|I_n coordinate|first three D4 coordinates|binary residual S|
|---|---|---|---|---|
|isolated|1,2,4|0|0|q|
|B_(n+1)+B3|2|e1|0|1|
|B_n+B4|2|0|(-1,0,0)|1|
|B_n+F4|4|0|(-1,-1,-1)|1|

Here the last D4 coordinate is `-y`, and

    y^2-phi*x^2=S,
    0<x, y>0,       phi*x/y<=1,
    x^2<=S,         x'^2<=phi*S.                  (3)

The two connected q2 profiles require the displayed product1 at a
long tail. The F4 extension has q4 and product2 with a3; q2/q4
cannot be interchanged arbitrarily. Other length assignments either
give a second4-bond in a connected finite diagram or a nonintegral
pairing. The old omitted-wall acuteness and several divisibility
conditions are harmlessly omitted from the upper bounds (3).

## 4. Tiny exact census and actual endpoint

Writing x=a+b phi, the two strips in(3) give

|S|complete positive x list|
|---|---|
|1|1|
|2|1, phi-1|

Indeed b=(x-x')/sqrt(5) can only be0 or1. For S1, x=phi-1 violates
the conjugate bound. For S2 it is admitted along with1; every other
integer a is outside one of the two bounds.

At S2 and x1, `y^2=phi+2` has rational norm5 and is nonsquare.
At S2 and x=phi-1, y=phi. At S1 and x1, y=phi again.

For either isolated q1 or q2, D4 requires y in2O, so these solutions
are excluded. For the B_(n+1) profile the same y-even condition fails.
For the B4 profile D4 instead requires `y=1 mod2`; phi is not1 in
O/(2), so that profile is excluded too.

At the F4 profile q4 forces x in2O. Write x=2u. Its residual S1
gives `0<u<=1/2`, `|u'|<=sqrt(phi)/2`, which contains no O element.
Thus this profile is empty.

It remains to consider isolated q4. Now x=2u, (3) gives `0<u<=1`
and `|u'|<=sqrt(phi)`; the only possibility is u1. This reconstructs
`x2,y2phi`, exactly r, which satisfies every omitted condition.
Therefore r is the unique first simple wall of the FULL arrangement.
Its actual endpoint has whole finite diagram `B_n+B3+A1(r)`.

## 5. The next gate is geometric, not a prefix loop

The r-perpendicular lattice is an index-ideal(2) submodule of an odd
unimodular diagonal lattice. Projection of `(a;u;w1,w2,w3,w4)` onto
r-perp has integral diagonal coordinates

    (t=phi*a+w4; u; w1,w2,w3),

with quadratic form `<-phi>+I_(n+3)`. This does NOT identify induced
arrangements. For example the actual old a4 mirror intersects r's
hyperplane, and its projected normal has norm

    a4^2-(a4.r)^2/r^2=4-phi^2=3-phi,

of rational norm5. Its ray is not a norm1/2 reflective ray of the
unimodular diagonal lattice. Thus the full-facet shortcut used for
the determinant-three architecture has a genuine obstruction here.
No statement that a4 is adjacent to r in the final chamber is needed
for this whole-arrangement observation.

A small repo/literature screen found no literal R_5 exclusion; that
is not a worldwide classification or novelty theorem. Before any
execution or further host work, the remaining gates are an exact
dyadic-module/known-host comparison and an actual outgoing-edge or
compactness argument. No unbounded prefix, compactness claim, global
period, new manifold, or dimension improvement is asserted.
