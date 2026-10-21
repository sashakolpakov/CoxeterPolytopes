# Mixed right-angled three-faces: a local ideal-edge corner obstruction

19 October2026. **New written candidate, awaiting independent review.**
No scientific job, graph/polyhedron census, optimizer, checker, Lean or
build ran. All existing manuscripts, checkers and certificates remain
unchanged. The unrestricted mixed right-angled dimension ceiling remains12.

The new local observation excludes the entire previously retained
F8/c3/p6 all-ideal-edge corner exception. It removes the finite--finite
edge terms from two current universal carrier-error estimates. A further
conditional corollary needs only ONE all-ideal-edge finite vertex rather
than global absence of finite--finite edges. None of these statements is
an unconditional dimension11 theorem.

## 1. A three-unit polygon-excess lemma

Let X be a finite-volume right-angled hyperbolic three-polytope, with
F facets, c ideal vertices and p finite vertices. Put nu=F+c. For a
polygonal facet A let u(A),v(A) be its finite and ideal vertex counts
and write

    e(A)=u(A)+2v(A)-5.

Hyperbolic area and integrality give e(A)>=0. Finite and cusp link
degrees are3 and4; Euler gives p+2c=2F-4. Therefore

    Delta(X)=sum_A e(A)=3p+8c-5F=F+2c-12.        (1)

Call a triangular facet with one finite and two ideal vertices bad.
Use the reviewed attachment lemma: if both facets on an ideal--ideal
edge are bad triangles, the whole X is the canonical P3. This is an
actual normal/convex-hull theorem, not a combinatorial assumption.

**Lemma.** If X is not P3 and a finite vertex q has three incident
edges ending at ideal vertices, then

    Delta(X)>=3, equivalently F+2c>=15.          (2)

Proof. Use the intrinsic Lorentz model of X, q^2=-1, with orthonormal
outward incident normals n1,n2,n3. Its three distinct ideal neighbors
are ui=q-ni. Let j be the number of bad triangular facets at q.
Every other incident facet contains q and TWO of these distinct ideal
vertices. Since it is not the triangle on those three points, it has
at least four vertices and e>=1. This gives3-j distinct positive-excess
facets through q.

For each bad triangle on directions i,k, its ideal--ideal edge has the
actual opposite wall

    m_ik=q-ni-nk.                               (3)

The facet on the other side of this edge contains its two ideal
endpoints. It cannot be another bad triangle, by the attachment lemma
and X!=P3. Its excess is consequently at least1: a polygon with two
ideal vertices has e=0 only when it is precisely a bad triangle.
For distinct pairs the walls(3) are distinct, by independence of the
incident orthonormal normals. Also q.m_ik=-1, so none of these partner
facets is incident to q. They give j further distinct positive-excess
facets. All3 facets just counted are distinct, proving(2).

This argument counts polygon excess, not the number of triangular
facets or a global ideal-edge graph. It allows other finite vertices,
finite--finite edges elsewhere, and arbitrarily many further facets.

## 2. The exact nu<=11 consequence

Suppose q has three ideal incident edges and nu<=11. If X=P3 there is
nothing to prove. Otherwise (2), together with p+4c=2nu-4 and p>=1,
gives

    nu+c>=15,
    4c<=2nu-5.

At nu<=10 these are incompatible. At nu11 they force c4,F7,p2.
Because q has no finite-edge neighbor and there is only one other finite
vertex, the other finite vertex also has only ideal-edge neighbors.
The reviewed seven-facet/four-cusp rigidity argument therefore identifies
X as the reflected double D3 of P3. In particular its finite corner
graph is an induced two-edge path, its missing-pair two-face is an
actual Q22, and the two finite vertices have cosh distance3.

For clarity, the previously allowed candidate F8,c3,p6 at nu11 fails
already in(1)--(2): its total polygon excess is2, while ANY all-ideal-edge
finite corner forces at least3. In particular it cannot support an
unbacked induced-path corner. No classification of its remaining
finite-edge corners is asserted or required.

Thus we obtain the unconditional local statement:

    At an all-ideal-edge finite corner of a nu<=11 three-face,
    that three-face is P3 or D3.                  (4)

The word local matters: global finite--finite edges were never excluded.

## 3. Stronger universal error bounds

Use the notation of
`right_angled_mixed_symmetric_quad_charge_research_20261018.md`.
At an ambient finite vertex, a path in its bad-triangle graph uses
three ideal directions. Its three-face has an all-ideal-edge corner.
By(4), every such path at nu<=11 is backed by the actual D3 symmetric
quad. Hence EVERY unbacked path has nu>=12.

Let E12=sum_(nu>=12)(nu-11), and let M(q) count unbacked induced paths
at q. For any non-P3 three-face let Z count its induced-path finite
corners. Each contributes three finite/ideal edges and two bad triangles
with distinct ideal/ideal edges. The attachment lemma prevents reuse of
these edges by another bad triangle. Thus

    E_FI>=3Z, E_II>=2Z, 4c=E_FI+2E_II>=7Z,
    p>=Z, p+4c=2nu-4,
    4Z<=nu-2.                                  (5)

For integral nu>=12, this gives Z<=2(nu-11): at nu12, Z<=2; at nu>=13
even the real bound (nu-2)/4<=2(nu-11) suffices. Summing yields the
new unrestricted estimate

    sum_q M(q)<=2E12.                           (6)

There is NO finite--finite edge error term. In particular the old
nu11 contribution (1/5)*binom(n-1,2)*e_FF can be removed.

For the companion count, let A3(X) count all-ideal-edge finite corners
of three-faces X other than P3 or D3. Statement(4) eliminates every
nu<=11 contribution. At nu>=12,

    3A3<=4c, p>=A3,
    A3<=(nu-2)/2<=5(nu-11).

Hence its unrestricted sum satisfies

    sum_(X neither P3 nor D3) A3(X)<=5E12.       (7)

The old companion term (2/3)*binom(n-1,2)*e_FF can also be removed.
These are universal inequalities at mixed corners. They do not assume
global absence of finite--finite edges or any nu upper bound.

## 4. One all-ideal finite corner already suffices under nu<=11

**Conditional corollary.** In dimensions9--12, a finite-volume
right-angled polytope with every three-face satisfying nu<=11 has NO
finite vertex all of whose incident edges end ideally. Equivalently,
every finite vertex then has at least one finite--finite incident edge.

Proof. Assume one such vertex q exists. Every three-direction face
there is P3 or D3 by(4). Consequently every induced triple of its
bad-triangle graph G(q) is a triangle or two-edge path. The complement
of G(q) is a matching: two incident missing edges would give a triple
with at most one bad edge. Let r be the number of missing pairs.

Use all actual P3 oppositions, with weight1, and one actual symmetric
quad move for each missing pair, with weight1/2. The existing cross-edge
identities are exact at both kinds of actual moves. Across a P3 move
on a triple S, each untouched direction has at least two bad neighbors
in S, since the complement is a matching. Its edge remains ideal;
the three new directions are ideal in the opposite P3 itself. Across
a missing-pair quad move, every direction has a bad cross-neighbor
outside/inside that pair (n>=9), again preserved by the exact quad
identities. Thus the opposite vertex has all ideal incident edges in
either case. By(4) its graph again has matching complement. This proves
that the finite move component starting at q is CLOSED within the
all-ideal-corner class; it does not assume that class globally.

On this closed finite component the existing exact displacement sum
has time/spatial coefficients

    a=binom(n,3)-r(n-3),
    b_un=binom(n-1,2)-r            (n-2r directions),
    b_pair=binom(n-2,2)-r+2        (2r directions).

Its causal margins a^2-sum b_i^2 are

    n9: 7r^2+6r;
    n10: 1440-168r+15r^2;
    n11: 4950-488r+25r^2;
    n12: 12100-1008r+37r^2.

They are nonnegative on integral0<=r<=floor(n/2), strictly positive
except n9,r0. In that endpoint the displacement is nonzero future
null; in every case a>0. All vertex displacement sums are therefore
nonzero future causal. Undirected weighted edge cancellation makes
their sum zero on the component, a contradiction. This is precisely
the existing causal calculation, now with a locally proved closed
domain in place of the old global no-finite-edge hypothesis.

It follows in particular that e_FF>=f/2 under global nu<=11 in these
dimensions. This is still a conditional structural constraint, not a
proof that the entire polytope is absent. Finite-edge carrier geometry
and the graphs with k<n ideal directions remain unresolved.

## 5. Current frontier and verification boundary

The reviewed ideal5/6 counts63/211 and275/6106 imply the mixed fractions
I5/N5<4/47 and I6/N6<15/248. They do not by themselves exclude a mixed
dimension12 system: zero entirely-ideal fractions are allowed. The
present result improves a different interface, eliminating the entire
low-nu carrier error without an ideal-face proportion assumption.

Even with(6), the full triangle/quad displacement need not be causal
for every graph: K8 disjoint union K2 disjoint union K2 remains an
abstract obstruction with no induced paths. No realizable example is
asserted. In the unrestricted case E12 can still be positive, and no
upper bound proved here makes its error smaller than a universal
causal margin. The dimension ceiling is therefore still12.

Written inherited inputs are the actual bad-triangle normal formula,
the attachment theorem, actual P3/D3 closure, the exact cross-edge
identities and the already reviewed matching-complement causal rows.
They are contained in the weighted/rigidity, actual-move-support and
symmetric-quad packets of18October. Their earlier sources remain
unchanged; the removed nu11 errors were conservative allowances, not
claimed existing geometries. No new finite-model calculation or complete
mixed-polytope census is substituted for these geometric proofs.
