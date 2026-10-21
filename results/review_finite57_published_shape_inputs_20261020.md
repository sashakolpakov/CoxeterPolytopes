# Primary-source audit of finite57's shape inputs

20 October 2026. Audited manuscript: `manuscripts/finite_volume57.tex`,
SHA-256 `c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0`,
lines 795–831, together with the curvature definition at 649–655.
Read the cited
[original English Khovanskii paper](https://www.math.utoronto.ca/askold/1986-Faa-1-english.pdf),
*Functional Analysis and Its Applications* **20** (1986), 41–50.
No canonical edit, scientific worker, Lean job or TeX job was run;
no PDF asset was saved. The separately prepared formal assembler was
not treated as evidence for this audit.

**Finding:** no gap in the published-input scope, incidence signs, or
internal derivation of `377*v_f + sum kappa + 4205` was found. The
following distinctions identify what the cited paper supplies and what
the manuscript still derives internally.

## Published statements actually checked

* Theorem 2, p.43 (PDF p.3), identifies the number of vertices having
  exactly `j` descending incident edges with coefficient `h_j` of the
  simple polytope's H-polynomial. Its corollaries give symmetry and
  positivity.
* Theorem 6, p.46 (PDF p.6), counts faces avoiding a level hyperplane by
  assigning below-section faces to their maxima and above-section faces
  to their minima.
* Section 5.1 and the **proof** of Theorem 9, p.49 (PDF p.9), define the
  projectivized vertex-cone figure `A(v)` and induced section `L(v)`.
  The proof explicitly contains both identities needed here:
  `F_(0,k)(A)=sum_v F_(k-1)(A(v))` and
  `2F_k(A)=sum_v F_(k-1)(A(v),L(v))`.
* Section 6, pp.49–50, explicitly treats finite-volume reflection
  chambers: they are simple at edges, with local cones over products
  of simplices, and may be nonsimple at ideal vertices. It also explains
  passage through the Klein model. Thus Theorem 9 is not restricted
  here to chambers with only finite vertices.

The manuscript uses identities from the reduction theorem's proof,
rather than mistaking its average-ratio conclusion for a new identity.

## Independent reconstruction of the manuscript application

Work in the compact Euclidean Klein closure `A`. Here `I` counts
vertex/two-face incidences, including ideal vertices, so
`I=sum_F s(F)`. A two-face through `v` gives exactly one edge of its
57-dimensional vertex figure. Consequently `I=sum_v E_v`.

Restrict the affine function `L-L(v)` to a bounded section of the
vertex cone. An edge of that figure is uncut exactly when the two
outgoing polygon-edge directions have the same sign under `L-L(v)`.
This means that `v` is an extremum of the corresponding convex polygon.
Each polygon has precisely two such vertices; therefore
`2f_2=sum_v E_v^0`. This also directly checks the factor two and the
sign convention without importing a ratio estimate.

Choose all small bounded figures first. The ambient function can avoid
the finitely many conditions of being constant on an original edge or
on an edge of any chosen figure. Its induced zero section avoids every
figure vertex because it is nonzero along every original edge. Thus
the genericity needed for all restrictions can hold simultaneously.
An empty induced section at an ambient extremum is allowed. No Lorentz
unit-normal convention or metric-face Coxeter assertion enters this
affine argument.

Orient a figure by increasing `L`. Below its zero section, assign an
uncut edge to its upper endpoint; above, assign it to its lower endpoint.
An index-`j` vertex receives respectively `j` or `57-j` edges. Hence
`E_v^0 >= sum_j h_j min(j,57-j)=G_v`.

A finite vertex figure is a 57-simplex. A cusp figure is projectively
equivalent to its bounded horospherical product of simplices, so is
simple even when the ambient cusp is nonsimple. Its H-polynomial is
the manuscript's `product_i(1+...+t^(r_i))`, and it has
`V=product_i(r_i+1)` vertices and `E=57V/2` edges. This supplies exactly
the interpretation needed for `kappa=(5G-57V)/2`.

At each of the distinct global minimum and maximum, `E_v^0=E_v`.
Counting edges by their upper endpoints gives

`E_v-G_v=sum_(j=29)^57 h_j(2j-57) >= 1+3+...+57 = 841`.

Every coefficient is at least one for these simplex/product polynomials.
Thus the two extrema add `(5/2)*(841+841)=4205` beyond the sum of
baseline curvatures, without counting any other vertex twice. For a
finite figure, `G=2*(1+...+28)=812` and `E=58*57/2=1653`, giving
`(5/2)G-E=377`. Finally

`Dfive=5f_2-I=sum_v((5/2)E_v^0-E_v)`

proves precisely the manuscript's claimed lower bound.

## Verification boundary

The published incidence/index identities and reflection-chamber local
structure were checked in the original source. The assignments,
simultaneous generic choice, simplex/product interpretation, numerical
constants and extrema correction above are the manuscript's internal
counting argument, independently reconstructed here. This audit does
not verify the separate payment rule, cusp cost envelope, finite-table
coverage or any new formal assembler execution.
