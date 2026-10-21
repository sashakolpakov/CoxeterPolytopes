# The sharp finite-edge triangle coefficient28/3

**Complete internally checked lemma.** Let P be a finite-sided,
finite-volume acute hyperbolic Coxeter polytope. An actual six-simplex
face having at least two finite vertices has at most four ideal vertices.
In particular, its two-finite/five-ideal case is impossible.

Together with the frozen higher-dimensional simplex bounds, this gives

\[
 \boxed{\quad2(T_0+T_1)\le {28\over3}E_{ff}
       ={14\over3}n v_f-{14\over3}E_{fi}.\quad}       \tag{1}
\]

The constant is attained by a classical finite-volume H9 Coxeter simplex
with seven finite and three ideal vertices. It is therefore sharp for
this universal finite-edge inequality. This is not an ambient dimension
bound or a new hyperbolic example. The demands from triangles with two
or three ideal vertices and quadrilaterals remain separate.

## Original stars and the539 cases

The frozen [simplex-star proof](fresh_multi_ideal_triangle_scope_20261017.md)
provides seven black facet roots B and a positive-definite supporting
white block W. The whole Gram is indefinite, while deleting any black
root gives the PSD Gram at the corresponding actual vertex; it is
positive definite exactly when that vertex is finite. All roots may be
positively rescaled to norm two. Induced face angles need not be Coxeter
angles, and ideal vertices of the ambient chamber need not be simple.

Choose an indefinite principal subdiagram L of minimum order in this
star. Every such diagram contains all seven blacks. The inherited
minimal-indefinite classification bounds its order by ten. The
[exact frozen census](fresh_high_rank_quasilanner_census_20261017.md)
contains three order-seven, four order-eight, four order-nine and three
order-ten types, giving

\[
 3+4\binom81+4\binom92+3\binom{10}3=539.
\]

The classical input remains the quasi-Lannér classification, recorded
with primary references in that census and
[Chein's original table, printed pp.15–16](https://www.numdam.org/item/M2AN_1969__3_3_3_0.pdf).
The current package does not rerun its large enumeration.

The [frontier certificate](fresh_six_simplex_ideal_frontier_20261017.md)
reconstructs all539 choices. In175 cases at least two black deletions
are connected affine: their positive kernels force every outside white
orthogonal to their union, namely all of L. In290 cases exactly one
black deletion is affine: the outside-white Schur correction changes
only that black diagonal. Its exact reciprocal-inverse thresholds and
all minimum ties are checked. These465 cases have ideal-count envelopes
`{2}`, `{3}`, `{1,2}`, `{1,3}`, or `{1,7}`. With at least two finite
vertices they already have at most three ideals.

It remains to treat74 cases with no affine BLACK deletion: four affine
rank-six black blocks and70 positive-definite rank-seven black blocks.

## Complete direct profiles and arbitrary indirect whites

Every remaining L has an affine deletion F=L-h with h WHITE. F is a
connected affine diagram. If a white u outside F touches F, the positive
affine kernel has strictly negative pairing with u, so F+u is indefinite.
Every proper subset of F+u has order less than L and is PSD by the
minimum-order choice. Therefore F+u is one of the complete same-order
quasi-Lannér types. This proves a finite direct-white alphabet without
restricting the number or rank of all supporting whites.

The frontier source enumerates every affine deletion in that catalogue
and every label-preserving graph isomorphism onto the literal F. Thus
all symmetrically distinct attachments are retained. Let U be the whites
outside F that touch it, including the mandatory root h. Let Z contain
the other whites, which are right to all of F. Z can have arbitrary
size, components and tails. Its Gram is positive definite because Z is
part of the supporting white block. Eliminating Z gives

\[
 D=R^tG_Z^{-1}R\succeq0,\qquad
 K=G_U-D\succ0.                                  \tag{2}
\]

Both D and the absolute values of the native offdiagonal products are
entrywise nonnegative: a positive-definite Gram with nonpositive
offdiagonals has a nonnegative inverse. In particular K has diagonal
at most2 and offdiagonal at most0. Eliminating the positive-definite
block F-b at each BLACK deletion b gives the exact necessary condition

\[
 K-C_b\succeq0,\qquad C_b=P_b^tG_{F-b}^{-1}P_b.    \tag{3}
\]

That vertex is ideal exactly when (3) is singular. The frontier records
all literal C_b. No assumption on the coefficient field of the indirect
whites is used in (2)–(3).

Since (3) is a PSD matrix with nonpositive offdiagonals, the entrywise
larger matrix `2I-C_b` is PSD too. Indeed for any real x,

`x^t(2I-C_b)x >= |x|^t(2I-C_b)|x| >= |x|^t(K-C_b)|x| >= 0`.

If a profile of maximal diagonal cost c occurs t times, the all-ones
vector on those copies gives `2t-t^2 c>=0`, hence `t<=floor(2/c)`.
All these costs are positive. The [theorem checker](../scripts/fresh_finite_edge_cusp_tradeoff_28over3_20261017.py)
enumerates every multiset within these proved multiplicity bounds and
tests every relaxed black-vertex matrix. This is not an assumed cutoff.
The resulting145 possibilities consist of111 single-profile models,
31 distinct-profile pairs, and three triples. Every rejected multiset
has a recorded indefinite black deletion.

## Single profiles and the651 five-vertex tests

For a single profile, repeated any permitted number of times,
`C_b=c_b 11^t`. Positive definiteness of K shows that `K-c_b 11^t` can
be singular only at the largest c_b allowed by all vertices. The checker
records every maximum and tie. There are at most four tied black
vertices, except one simultaneous all-seven alternative. Thus exactly
five ideal vertices is impossible in all111 models.

For each of the31 distinct-profile pairs write

\[
 K=\begin{pmatrix}x&z\\z&y\end{pmatrix},\qquad t=\det K.
\]

Its black-vertex determinant is

\[
 \det(K-C_b)=t-x(C_b)_{22}-y(C_b)_{11}
                  +2z(C_b)_{12}+\det C_b.          \tag{4}
\]

All coefficients here are rational, as checked from the literal
matrices. The21 choices of five black vertices give651 exact linear
systems in x,y,z,t. Their complete outcomes are:

| Reason | Systems |
|---|---:|
| The determinant equations are inconsistent |628|
| They force at least six black determinants zero |21|
| The remaining determinant identity forces z>0 |2|

The second outcome leaves at most one finite vertex and is excluded by
the theorem's hypothesis. The last two systems are the two symmetric
affine-A6 configurations. Their linear equations force

\[
 x=y={7t\over12},\qquad z={7t\over12}-{3\over2}.
\]

Substituting `t=xy-z^2` gives t=3 and z=1/4, contradicting K's
nonpositive offdiagonal. The saved certificate contains every five-set,
coefficient matrix, exact ranks and forced-singularity closure; these
are exact linear-algebra exclusions, not numerical feasibility tests.

## The three triples

Write `G_U=2I-N`, where N is entrywise nonnegative with zero diagonal.
For a nonnegative null vector v of `H_b=2I-C_b`, PSD of
`H_b-N-D` implies `(N+D)v=0`. Since every entry of N+D is nonnegative,
every row and column meeting the positive support of v must vanish.

The affine-E6 triple has null vectors `(1,1,0)`, `(1,0,1)` and `(0,1,1)`
at three black vertices. Their supports cover all coordinates, forcing
N=D=0. Exactly those three vertices are ideal.

The two remaining triples are L22 with internal whites `{4,8}` or
`{7,8}`. They contain two equal profiles and one other profile. A
nonnegative null vector supported on the repeated pair forces N=0 and
all rows of D except the singleton coordinate to vanish. Hence
`D=tau e e^t`, tau>=0. The complete black-vertex threshold lists are

| Internal whites | Black labels | Thresholds |
|---|---|---|
| `{4,8}` | `0,1,2,3,5,6,7` | `5/4,1,4/3,3/2,7/6,1,1/2` |
| `{7,8}` | `0,1,2,3,4,5,6` | `5/4,1,7/6,1,1/2,4/3,3/2` |

Each has unique minimum1/2, distinct from its pre-existing ideal vertex.
Thus each has at most two ideals. At an already singular matrix the
threshold is computed by grounding a coordinate where its unchanged
kernel is nonzero; the corrected coordinate is orthogonal to that
kernel, so the grounded inverse gives the same scalar PSD condition.
This treats the entire half-line of possible indirect corrections.

The465 easy cases and145 models for the74 hard cases exclude exactly
five ideals when at least two vertices are finite, completing the
six-simplex assertion with arbitrary supporting whites.

## Coefficient and sharpness

The [finite-ended simplex identity](fresh_finite_ended_simplex_ports_20261017.md)
assigns each finite edge to a unique maximal actual simplex block with
f>=2 finite vertices and i ideal vertices. Its twice-triangle payment
per finite edge is `2i+2(f-2)/3`. The exact remaining maxima are:

| Simplex dimension | Maximum |
|---|---:|
| At most5 |8|
|6, with i<=4 |26/3|
|7, with i<=3 |8|
|8, with i<=3 |26/3|
|9, with (f,i)=(7,3),(8,2),(9,1) |28/3|

The higher-dimensional domains are inherited from the
[coefficient10 theorem](fresh_finite_edge_cusp_tradeoff10_20261017.md),
whose saved replay was rechecked at the start of this work. Summing
the block inequalities proves (1). The corresponding residual budget is

\[
 2(T_2+T_3)+Q\ge(b_n-\tfrac{14}{3}n)v_f
        +\tfrac{14}{3}E_{fi}+\alpha_n\sum_v E(L_v).
\]

Sharpness is literal, rather than a formal maximizer. Classical rank-ten
diagram L24 has three affine deletions and seven finite deletions. The
checker reconstructs its Lorentz Gram, inverse, strict timelike center,
all ten future causal vertex rays, and the positive affine kernels of
its three full-rank cusp links. It therefore defines a finite-volume H9
Coxeter simplex with seven finite and three ideal vertices. It has
`T0=35`, `T1=63`, and `E_ff=21`, giving `2(35+63)/21=28/3`.
This is an already classified simplex, not a new construction claim.

## Verification and limits

The [exact output](fresh_finite_edge_cusp_tradeoff_28over3_20261017.json)
binds the frontier, frozen census, arithmetic engine and preceding
geometric proofs. Generation and fresh saved replay reconstruct every
small calculation. The [independent internal review](review_finite_edge_28over3_20261017.md)
passed both mathematical/source review and a fresh exact generation and
saved replay. Its separate Fraction Gaussian engine checks all111
single-profile envelopes,651 determinant systems,three triple cases,
and the sharpness simplex without importing the producer arithmetic or
SymPy. It treats the539-case catalogue and profile isomorphisms as
reviewed inputs rather than independently repeating that enumeration.
A second conceptual review confirmed the minimum-order completeness,
arbitrary-white Schur reduction and nonnegative-kernel arguments.
The frozen producer JSON retains its historical candidate status string;
the completed review establishes the current status. The
white-profile classification is a finite consequence of minimum order;
the arbitrary-Z Schur argument is the all-size proof that prevents a
hidden finite-tail cutoff. No new ambient dimension conclusion follows
until the remaining multi-ideal triangles and quadrilateral demand are
controlled. This is internal verification, not external refereeing or
proof-assistant formalization. **FINAL/FROZEN.**
