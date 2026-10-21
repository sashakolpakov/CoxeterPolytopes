# Finite-vertex parity, a five-dimensional seed, and inherited face covers

20 October 2026. New paper lemma for the mixed-dimensional research lane;
no checker, optimizer, census, Lean or build was run. The large-dimensional
cover consequences below explicitly use the separately proved ideal
dimension-four theorem in
`results/ideal_right_angled_dimension4_proof_20261020.md`. They do not
establish a new mixed dimension bound.

## 1. A local parity covering lemma

**Lemma.** A right-angled finite-volume j-polytope having a finite vertex
has at least j finite vertices when j>=5. In particular, a mixed
five-polytope has at least five finite vertices.

At a finite vertex q the j incident facets are distinct and mutually
orthogonal. Every three-subset S of their omitted labels determines an
actual three-face Y_S through q. Every right-angled three-face has an even
number p of finite vertices: Euler and the finite/ideal degrees give

    p+2s=2F-4.

Consequently Y_S must contain a finite vertex other than q.

For any other finite vertex r define its nonempty missing-label set T(r):
the incident q-facets which do not contain r. It is nonempty because all j
facets intersect only at q. The face relation is exact:

    r belongs to Y_S if and only if T(r) is a subset of S.       (1)

Choose one label i(r) from each T(r). If there were at most j-3 other
finite vertices, at least three labels would be unchosen. Their
three-subset S contains no i(r), so (1) shows that Y_S has only the one
finite vertex q, contradicting parity. Thus there are at least j-2 other
finite vertices, giving the preliminary lower bound j-1.

Suppose equality j-1 held. There are exactly j-2 other finite vertices.
Their chosen labels must all be distinct, since otherwise the same
three-unchosen-label argument applies. Rename the vertices r_i and labels
i=1,...,j-2; call the two unchosen labels a,b.

The face Y_{i,a,b} can contain no r_l with l!=i, because its chosen label
l is absent. Parity forces it to contain r_i. Thus

    T(r_i) is a subset of {i,a,b}.                             (2)

For distinct i,l, the face Y_{i,l,a} can contain only r_i,r_l besides q.
Its finite count must be even, so it contains exactly one of them. By
(2), if beta_i denotes the bit `T(r_i) subset {i,a}`, this says

    beta_i+beta_l=1                                           (3)

for every distinct pair. Since j-2>=3, three such bits cannot satisfy all
three pair conditions. This contradiction excludes equality j-1 and
proves the lemma. No classification of the three-face or assumption on
the types of its polygons was used.

## 2. Odd-dimensional parity and the seed six

In dimension j>=3 a finite vertex has degree j, whereas an ideal vertex
has degree 2^(j-1), the number of vertices of its cubical link. Edge
handshaking therefore gives

    j*f_finite+2^(j-1)*f_ideal=2*f_edges.

For odd j, the finite-vertex count is even. In particular the previous
five-dimensional lower bound five improves to **six**. This is also a
consequence of odd-dimensional truncation Dehn--Sommerville symmetry,
but no polynomial computation is needed for this parity statement.

## 3. Consequence of excluding entirely ideal dimensions above four

Every actual face of dimension at least five is finite-volume,
right-angled and has a finite vertex: otherwise it would contradict the
separate ideal dimension-four theorem. Its finite vertices are simple.
Let v_j be a valid lower bound on their number and M_j a valid facet
minimum in dimension j. Every facet of a j-face, j>=6, has at least
v_(j-1) finite vertices. Counting finite-vertex/facet incidences gives

    j*f_finite >= M_j*v_(j-1),
    v_j=ceil(M_j*v_(j-1)/j)                                  (4)

as a valid recursive lower bound, with further rounding up to an even
integer in odd dimensions. The seed from Sections 1–2 is v5=6.
Using the inherited facet bounds

    M6=27, M7=46, M8=82, M9=218, M10=1923

recorded in the earlier mixed-incidence packet, (4) yields

|j|5|6|7|8|9|10|
|---|---:|---:|---:|---:|---:|---:|
|v_j|6|27|178|1825|44206|8500814|

These are necessary finite-vertex counts, not proposed examples. They
improve the seed obtained just by doubling a polytope with one finite
vertex; the parity proof avoids relying on doubling for this step.

## 4. The exact dimension-ten aggregate rows

In an ambient ten-polytope let f be its total finite-vertex count and a_j
the number of j-faces. A finite vertex belongs to exactly binom(10,j)
j-faces, so the new cover inequalities are

    6*a5<=252*f,
    27*a6<=210*f,
    178*a7<=120*f,
    1825*a8<=45*f,
    44206*a9<=10*f,
    8500814<=f.                                              (5)

The last row is the j=10 version of (4), not an incidence with an
unaccounted family of faces. In a normalization f=1 it becomes
`8500814*u<=1` with u=1/(the unnormalized finite-vertex count).

All rows require the actual nonexistence of entirely ideal five-faces,
not a statement that mixed polytopes inherit only mixed facets without
that theorem. They are suitable additional inputs to the joint-log
relaxation after independent review; no relaxation was executed here.
