# Twenty-eight vertices exclude an entire regular pentagon

20 October2026. Complete written candidate, awaiting independent review.
No enumeration, scientific execution, new worker, manuscript or formal
verification. The published count of two C28 fullerenes is an explicit
classification input. No geometric existence is inferred from a formal
Gram matrix; all metric arguments start from an assumed actual compact
right-angled polyhedron and its metric uniqueness.

**Theorem.** A compact right-angled hyperbolic three-polyhedron with28
vertices has no entire regular right-angled pentagonal face.

Together with the already reviewed26-vertex result and L6 metric, this
would make every positive facet adjacent to D in compact RA4 have excess
at least five. It gives f>=65 and Vol>=49pi^2/6, not the120-cell minimum.
These consequences are candidates until this proof is independently
reviewed; the frozen C26 packet and global-charge source are unchanged.

Write phi=(1+sqrt5)/2, u=phi^-1 and L=arcosh(phi). Scalar products below
are Lorentz bilinear products; all facet normals are outward and unit.

## 1. The three exhaustive combinatorial cases

Trivalence and Euler give42 edges and16 faces. The reviewed polygon
capacity F>=2k+2 bounds every face size by seven. If a heptagon occurs,
equality holds. The FULL equality-pattern proof after Lemma `polygon`
in the current compact RA4 manuscript gives the Loebell L7 face lattice:
the unique remaining vertex of the complementary dual disk must join
the whole opposite7-cycle, with no chord allowed by the ear argument.
This is stronger than merely knowing that one extra face is disjoint.
It gives two7-gonal bases and fourteen lateral pentagons.

Compact3D metric uniqueness realizes the base rotation isometrically,
so both bases are regular. Their edge length L7 satisfies

    cosh L7=1+2cos(2pi/7)>1+2cos(2pi/5)=phi.

Every lateral pentagon has a base edge of this length. Thus no pentagon
of L7 is regular of length L. This resolves every case with a heptagon.

Otherwise there are exactly twelve pentagons and four hexagons. The
one-skeleton is a fullerene on28 vertices. The only imported enumeration
is that there are TWO such fullerene graphs, from Brinkmann, Goedgebeur
and McKay, *The Generation of Fullerenes*, arXiv:1207.7010v2,
Table2, PDF page15, row(nv,nf)=(28,16), column“fullerenes”=2:

https://arxiv.org/pdf/1207.7010v2 .

The next two sections give two explicit nonisomorphic dual spheres with
this degree data; hence the published count makes them exhaustive.
Their geometries are derived here, not imported from molecular lengths.

For each explicit graph automorphism, compact3D Andreev uniqueness
provides an isometry of an ASSUMED actual Q with that combinatorics.
These isometries form a finite group fixing an interior point: normalize
the Lorentz sum of the orbit of any interior point. In particular the
metric equations below concern actual facets and actual edge endpoints.
They do not constitute a compactness or incidence test for a partial Gram.

## 2. The tetrahedral graph: one shorter edge on every pentagon

Let i,j,k,l range over four indices0,1,2,3. The16 vertices of the dual
are H_i and P_ij, i!=j. Start with the truncated tetrahedron on its12
directed-edge vertices P_ij. Its four triangular faces are

    {P_ij: j!=i},          i fixed.

For each k the other three indices i,j,l give a hexagonal boundary

    P_ij,P_ji,P_jl,P_lj,P_li,P_il.

Insert H_k in this hexagon and join it to every boundary vertex. Retain
the four triangular faces. This is visibly a triangulated sphere with
28 triangles. Every P_ij has five neighbors: P_ji, the two other P_i*,
and the two H_k with k not i or j. Every H has six neighbors. There are
no H--H edges. Every permutation of the four indices is an automorphism;
it acts transitively on all twelve pentagons P_ij.

At the symmetry-fixed interior point, the spatial components of the four
H normals form a regular tetrahedron: their lengths and pair products
are equal, and their4-by4 Euclidean Gram has rank at most three. Distinct
H faces exclude the rank-one alternative, forcing pair products -1/3
after unit normalization. Thus for tetrahedral unit vectors v_i,

    v_i.v_j=-1/3 (i!=j),     H_i=(h,-R*v_i),
    P_ij=(a,A*v_i+B*v_j),    R^2=1+h^2.

The stabilizer of ordered(i,j) fixes its spatial plane, justifying the
two-coefficient expression for P_ij; transitivity makes a,A,B common.
Here h,a>0. Orthogonality to H_k, k not i,j, gives

    A+B=3ah/R>0.                                      (1)

The adjacent pairs P_ij,P_ik and P_ij,P_ji give respectively

    a^2=A^2-(2AB+B^2)/3,
    a^2=2AB-(A^2+B^2)/3.                              (2)

The unit norm of P_ij together with the first equality gives B^2=3/4.
Subtracting (2) gives A(A-2B)=0. The branch A=0 would give a^2<0, so
A=2B. Equation(1) selects B>0. It follows that

    a^2=7/4,      h^2=3/4,      R^2=7/4.              (3)

The actual edge P_ij intersect P_ji has endpoint third normals H_k,H_l
for the other two indices. For a compact RA3 edge, the product of its
two endpoint third normals is minus the cosh of its length. Therefore

    cosh(length(P_ij intersect P_ji))=h^2+R^2/3=4/3.

This is strictly less than phi. Every pentagon contains one such edge,
so none is regular. No statement about absence of fivefold symmetry was
used to exclude an accidental regular pentagon.

## 3. The second graph and its actual four-element symmetry

The second16-vertex dual sphere has hexagonal labels H0,H1,H2,H3 and
pentagonal labels A,B,C,D,E,F,G,J,K,M,N,P. Its28 triangles are

    H0 H1 A;  H0 H1 B;  H0 B C;   H0 C D;
    H0 D E;   H0 E A;   A E F;    A F G;
    A G H1;   H1 G J;   H1 J K;   H1 K B;
    B K M;    B M C;    C M N;    C N D;
    D N H2;   D H2 E;   E H2 F;   H2 F P;
    F G P;    G J P;    J H3 P;   J K H3;
    K M H3;   M N H3;   N H2 H3;  H2 H3 P.

For completeness, its cyclic vertex links are

    H0: H1,B,C,D,E,A;          H1: H0,B,K,J,G,A;
    H2: D,E,F,P,H3,N;          H3: J,K,M,N,H2,P;
    A: H0,E,F,G,H1;           B: H0,C,M,K,H1;
    C: H0,D,N,M,B;            D: H0,E,H2,N,C;
    E: H0,A,F,H2,D;           F: A,G,P,H2,E;
    G: A,H1,J,P,F;            J: H1,K,H3,P,G;
    K: H1,B,M,H3,J;           M: B,C,N,H3,K;
    N: C,D,H2,H3,M;           P: F,G,J,H3,H2.

Every edge has two triangle incidences; connectedness and these links
give a closed connected2-manifold, with Euler16-42+28=2. Thus it is a
sphere. Its H-induced graph consists of two edges H0H1,H2H3, distinguishing
it from the tetrahedral graph. This and the primary count-two input prove
exhaustiveness without a new fullerene search.

The permutations

    alpha=(H0 H1)(H2 H3)(A B)(C G)(D J)(E K)(F M)(N P),
    beta =(H0 H2)(H1 H3)(A N)(B P)(C F)(D E)(G M)(J K)

preserve all triangles, commute and have order two. Their pentagon
orbits are

    {A,B,N,P},      {D,E,J,K},      {C,F,G,M}.          (4)

By the actual-Q uniqueness bridge these are actual isometries. The proof
does NOT need to assume they are rotations rather than reflections, or
to classify their representation on the ambient spatial three-space.

## 4. Two orbits eliminated by exact pentagon propagation

We use a small strengthening of the elementary regular-pentagon lemma:
ANY TWO edges of length L in a compact right-angled pentagon force the
whole pentagon to be regular. Adjacent positions give the familiar
two-edge-chain uniqueness. For nonadjacent positions s0=s2=L, the
right-pentagon identity

    cosh(s0)=sinh(s2)*sinh(s3)

gives sinh(s3)=phi/sqrt(phi)=sqrt(phi)=sinh L. The adjacent pair s2,s3
then gives regularity. Reversing the cyclic order covers the other
nonadjacent choice. The identity and adjacent-chain proof are already
derived in the compact RA4 manuscript.

If any face in the first orbit(4) is regular, all four are regular.
Then C has two L-edges from its neighbors B,N, so C and its whole third
orbit are regular. If instead the third orbit is initially regular,
A has two L-edges from F,G, reducing to the preceding case. With both
of those orbits regular, D has two L-edges from C,N, so the middle orbit
is regular as well. All twelve pentagons would now be regular. H0 has
five consecutive L-edges, since its only hexagonal neighbor is H1.
This contradicts the three-consecutive-L-edge closure lemma for a
nonpentagonal compact RA polygon. Thus neither the first nor the third
orbit can contain a regular pentagon.

## 5. The remaining orbit forces duplicate actual facet normals

Suppose a face of the middle orbit is regular. All D,E,J,K are then
regular. In particular their actual unit normals D,E are orthogonal,
both orthogonal to H0,H2, and regularity gives

    <H0,H2>=-phi.

The four vectors D,E,H0,H2 are a basis: the first two give a positive
orthonormal plane and the latter two a Lorentz plane of determinant
1-phi^2=-phi. In the wall of a regular pentagon, the scalar product
of any two nonadjacent side normals is -cosh L=-phi. Using the cyclic
links of D and E above, respectively, determines the actual normals
C and A uniquely in this basis:

    A=-phi*D+phi*H0+H2,
    C=-phi*E+phi*H0+H2.                               (5)

For example A has products0 with E,H0 and -phi with D,H2; C has
products0 with D,H0 and -phi with E,H2. These constraints give(5)
directly, and the resulting vectors have norm one.

Put <H0,H3>=<H1,H2>=-k. Equality follows by alpha; k>=1 by the compact
actual-wall gap, because these hexagons are not adjacent. Also
<H0,H1>=<H2,H3>=0 and, by alpha, <H1,H3>=-phi. Put e=<H1,E>.
The actual orthogonality <A,H1>=0 and(5) give

    H1=-(k/phi)D+eE+kH0+(k/phi)H2.

Apply beta, which interchanges D,E and H0,H2, to get H3. Unit norm and
the product <H1,H3>=-phi now give the two exact equations

    e^2=1+phi^-3*k^2,
    e=(phi-k^2)/(2phi^-1*k).                           (6)

Eliminating e, using only phi^2=phi+1, gives

    (33-20phi)k^4-(8-2phi)k^2+phi^2=0.

Its two roots in k^2 are

    phi^4,          (6+7phi)/29.

The second is less than1 since phi<2. Thus k>=1 forces k=phi^2, and
(6) gives e=-phi. The Gram matrix of H0,H1,H2,H3 is consequently

    [ 1,     0,    -phi,   -phi^2 ]
    [ 0,     1,    -phi^2, -phi   ]
    [-phi, -phi^2,   1,       0   ]
    [-phi^2,-phi,    0,       1   ].

These H vectors span a nondegenerate Lorentz three-space U. Explicitly,
H0+H2 is timelike; the Gram determinant of H0,H1,H2 is
1-phi^2-phi^4<0, so they are independent. The displayed Gram makes
H0-H1-H2+H3 orthogonal to all H vectors with norm zero. A vector
orthogonal to the timelike H0+H2 has nonnegative norm, vanishing only
for zero, so that relation is zero. Hence dim U=3 and U-perpendicular
is a positive one-dimensional space.

Using(5),(6), both A and C have products with H0,H1,H2,H3 equal to

    (0,0,-phi,-phi).

They therefore have the same projection onto U. They are distinct,
since A-C=-phi(D-E) is nonzero. On a positive one-dimensional
orthogonal complement, a prescribed projection and unit norm allow
at most two vectors, and A,C are those two possibilities.

But G=alpha(C) is a THIRD actual unit normal with exactly the same
products: alpha swaps H0,H1 and H2,H3, preserving the displayed row.
Thus G must equal A or C. Both contradict the distinct actual facets
in the explicit graph. This excludes the middle orbit and proves
the theorem for the second fullerene.

## 6. Whole-polytope consequence, conditional on acceptance of this proof

An actual positive RA4 facet adjacent to D contains their entire
regular pentagonal ridge. Excess two is excluded by the L6 base-edge
metric; excess three by the frozen C26 theorem; excess four by the
theorem above. Thus every such boundary-positive facet has d>=5.
There is always a D facet by the reviewed farthest-facet theorem. A
non-all-D compact polytope has an actual D--positive adjacency by
connectedness, hence a facet with V>=30. Facet capacity gives f>=65.
Together with W>=17 in the nonzero case this gives

    Vol(P)>=pi^2*(65+17+16)/12=49pi^2/6.

The all-D120-cell satisfies these bounds too. No boundary facet of
excess five has been recognized here, and no universal120-cell
minimum or second-volume result follows from this argument.

## Unchanged inputs and provenance boundary

* `manuscripts/compact_ra4_volume_constraints.tex`, SHA256
  `afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f`:
  full polygon equality-pattern proof; actual-wall gap; compact3D
  metric uniqueness; pentagon identities and two-edge uniqueness;
  three-L polygon closure; facet capacity, W gap and volume formula.
* `results/compact_ra3_twenty_six_vertices_no_regular_pentagon_20261020.md`,
  SHA256 `1bbe6adc630ff2f44424f6ef3c573b3c71dcef4f25496def9b13fd35b2e4a03f`:
  used only for the excess-three case in Section6.
* `results/compact_ra4_farthest_dodecahedron_wall_radius_20261020.md`,
  SHA256 `d082e4a918d4504a59ac200f32bbb141f1f2cb52d701006eaea7ba5b0e5bcc22`:
  existence of a D facet, used only for Section6's RA4 consequence.

The sole external finite-classification premise is the cited primary
count-two for C28. Both graphs, their actual symmetries and the metric
obstructions are explicit above. No classification of all compact
RA4 polytopes, geometric completeness search, executable check or
formal geometry is claimed. All previously frozen sources are intact.
