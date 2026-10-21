# A nu>=12 three-face is forced in mixed dimensions eleven and twelve

19 October2026. **New written candidate; independent review pending.**
No computation, census, optimizer, Lean, build or source modification of
earlier evidence was used. This is an unconditional carrier-existence
obstruction in ambient dimensions11/12, not exclusion of those dimensions.
The unrestricted dimension ceiling remains12.

**Theorem.** Every finite-volume right-angled hyperbolic polytope in
dimension11 or12 contains an actual three-face X with

    nu(X)=facets(X)+cusps(X)>=12.

The earlier universal theorem forced only nu>=11. No hypothesis about
finite--finite edges is used here. The entire previous conservative
nu11 unbacked-path class is avoided by the new local excess lemma.

## 1. Local input and maximal ideal degree

Use the just-reviewed local theorem
`right_angled_mixed_nu11_all_ideal_corner_exclusion_research_20261019.md`,
SHA256
`fc047b78e9a8fc9eaf283bd5d689b69d72a8bbf0976a4765600ffbe8625955b0`:
at an all-ideal-edge finite corner of a three-face with nu<=11, that
three-face is P3 or its reflected double D3. Its proof charges three
distinct positive-excess polygonal facets and uses the actual triangle
attachment theorem. It does not require absence of finite edges elsewhere.

Suppose, toward contradiction, that every three-face of the ambient
polytope has nu<=11. For a finite vertex q let k(q) be the number of
its incident edges whose other endpoint is ideal, and put

    K=max_q k(q).

At q the bad-triangle graph has edges only between these ideal
directions. By the local theorem, every triple within the ideal set
is P3 or D3. Thus the graph induced on that set is a complete graph
minus a matching: every triple has two or three bad edges. The remaining
directions are isolated in the bad graph because their edges are finite.

We prove that **K<=8** in ambient dimensions9--12 under this global nu
assumption. Suppose instead K>=9. Start at a K-maximal vertex and use
only actual P3 oppositions on ideal-direction triples and the actual
symmetric quad on every missing ideal pair. Each missing pair has
an ideal third direction; its three-face is D3, so the quad really
exists. Different carriers of that pair still supply only one move.

## 2. The maximal-K class is closed under these actual moves

At a finite vertex use q^2=-1 and orthonormal outward normals ni. For
a P3 triple S the actual opposite is q'=2q-sum_(i in S)ni. Cross
bad edges between S and its complement are preserved exactly. Since
the ideal-direction graph has matching complement, each untouched ideal
direction has at least two bad neighbors in S. It remains ideal at q'.
The three new directions are ideal in the opposite P3 itself. Thus all
K old ideal directions persist under the transported direction labels.

For a missing pair i,j, every other ideal direction is a bad neighbor
of both i and j. The actual quad move replaces their normals by

    Ri=2q-ni-2nj, Rj=2q-2ni-nj,
    q'=3q-2ni-2nj.

Its exact cross identities are
q'-Ri-nl=q-ni-nl and q'-Rj-nl=q-nj-nl. Because K>=9, the other
ideal directions are nonempty; they make both new directions ideal,
and each of those untouched ideal directions remains ideal as well.
Again at least K ideal directions persist.

By maximality, neither move can create an ADDITIONAL ideal direction:
its opposite must have exactly K. Its ideal-direction graph again has
matching complement by the same local three-face theorem. Induction
therefore gives a finite move component consisting entirely of K-maximal
vertices. All moves used here stay within it, and reversal gives the
same set and weight. The directions outside its K ideal directions can
have finite edges; they are not assumed absent or unchanged globally.

## 3. Effective-dimension causal contradiction

Give each actual P3 move weight1 and each missing-pair quad move weight1/2.
At a vertex of the component let r be the number of missing pairs among
its K ideal directions. Summing displacements q'-q has the form

    D(q)=a*q-sum_(ideal i) bi*ni,
    a=binom(K,3)-r(K-3),
    b_un=binom(K-1,2)-r              (K-2r directions),
    b_pair=binom(K-2,2)-r+2          (2r directions).

These coefficients may also be obtained directly: removing r matching
edges erases r(K-2) triples; the r half-weighted quad moves restore r
to the time coefficient and one to each paired spatial coefficient.
The remaining ambient directions have coefficient ZERO. Accordingly
the Lorentz norm calculation is in effective dimension K; it does not
require the ambient dimension to equal K.

For K9,10,11,12 the margins a^2-sum bi^2 are respectively

    7r^2+6r,
    1440-168r+15r^2,
    4950-488r+25r^2,
    12100-1008r+37r^2.

For integer0<=r<=floor(K/2) these are nonnegative, strictly positive
except K9,r0. The last three quadratic functions decrease on the allowed
real intervals and have respective positive integer-domain minima
975,3135,7384. In all cases a>0 and D(q).q=-a<0; D(q) is a nonzero
future causal vector. This includes the K9,r0 null case.

But the sum of D(q) over the closed finite component is zero by
undirected weighted edge cancellation. A sum of nonzero future causal
vectors cannot vanish. This contradiction proves K<=8.

## 4. Existing universal triangle density finishes the carrier theorem

Let B count actual P3 three-faces, and f the finite vertex count. The
actual triangle correspondence gives

    2B=sum_q t(G(q)).

If K<=8, every graph triangle uses at most eight ideal directions, so
2B<=56f. The existing reviewed weighted-three-face theorem already gives
the weaker bounds sufficient here:

    ambient11: 2B>(7545/118)*f>56f,
    ambient12: 2B>(5030/59)*f>56f.

These are explicit consequences in the original weighted packet and
primary corner manuscript; the later stronger polygon-coupling lower
bounds are not needed. In these dimensions a finite vertex exists,
as also follows from the standard ideal right-angled dimension ceiling
and from those weighted ratio statements. Thus f>0, and the bounds
contradict2B<=56f. The supposition that all three-faces have nu<=11 is
false, proving the theorem.

The proof does not turn a maximum-energy vertex into a global average.
It uses the ordinary MAXIMUM NUMBER of ideal incident edges to obtain
a closed move class, then an exact vector cancellation and an existing
global count. The number K cannot increase along a move only because
it was chosen globally maximal; this is stated before invoking closure.

This supplies a genuinely stronger mandatory carrier, and focuses the
remaining unrestricted incidence problem on nu>=12 three-faces. It
does not quantify their density or prevent them from providing all
needed errors. The coefficient-one/two error estimates in the separate
nu12-sharpening note remain necessary constraints, not a completed
dimension12 contradiction or a new example.
