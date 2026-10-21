# Golden determinant-three scope and a uniform diagonal barrier

20 October 2026. Paper proof and literature scope only. No worker,
root search, Lean or TeX job was run. No new example or dimension bound.

## Exact published-input boundary

The own complete-edge/full-facet proof is
`results/golden_determinant3_g2_first_edge_and_full_facet_20261020.md`,
SHA256 `24870d50f72e768ea3dc4734680c49689b462da56af0850d6bb58e05e9dfa51e`.
For `L_n=<-phi>+I_n+A2`, n>=7, it identifies an actual FULL facet as
`M_(n+1)=<-phi>+I_(n+1)`, not merely a lattice of that determinant.

Kellerhals, *Hyperbolic Orbifolds of Minimal Volume*, Comput. Methods
Funct. Theory14(2014),465--481, §4.1.1, printedp474, equation(4.3),
explicitly identifies the literal form `M_j` and attributes reflectivity
only through j=7 to Bugaenko1984 [4]. The following paragraph identifies
the DIFFERENT dimension8 form of discriminant `-(1+sqrt(5))` from
Bugaenko1992, Theorem3.4. Author-hosted readable source:
https://homeweb.unifr.ch/kellerha/pub/CMFT.pdf
DOI: https://doi.org/10.1007/s40315-014-0074-y

This is an explicit later published attribution, **not** a claim that
the original1984 proof has been read. Its original MathNet article is
https://www.mathnet.ru/eng/vmumm3408 . The905KB PDF remained inaccessible
to the bounded web-reading attempts; no further access loop was made.
The weaker isolated phrase in Vinberg1985 is not used alone to infer
nonreflectivity.

Conditional on the exact literal cutoff just stated, the new L_n
family has no compact full reflection chamber for every n>=7. Indeed
the full coordinate facets of M_j are M_(j-1): intrinsic roots of
norm1/2 extend, and the only nonorthogonal crossing is a norm2 root
with coordinate product+-1, projecting to an integral norm1 root.
Thus compact M_j for j>=8 would imply compact M_8. This conclusion is
explicitly dependent on the published cutoff input, not an independent
new proof of that cutoff or a commensurability-class exclusion.

## A uniform reason not to try more positive diagonal coefficients

Let M be an admissible unimodular O_phi lattice, and let d be ANY
nonunit totally-positive element. In `R=M+<d>` the coordinate normal c
has norm d, is reflective, and its FULL induced arrangement is the
full integral reflection arrangement of M.

To prove this, write a primitive ambient normal as `w=v+k c`, v in M.
Its norm ideal `(q)` divides `(2d)`. A distinct mirror crossing c has
`v^2=q-dk^2>0` at BOTH embeddings: at the distinguished place this
is the crossing condition, and at the other it is positivity of M.
If `(q)` is a proper divisor of `(2d)`, every prime ideal of O_phi has
norm at least4, so `N(q/d)<=1`. A nonzero k would then give
`N(k)^2<N(q/d)<=1`, impossible.

In the remaining case unit normalization gives `q=2d`, because every
totally-positive unit is a square. A nonzero transverse k must satisfy
`|k|,|k'|<sqrt(2)`, which leaves only k=+-1 in O_phi. (The other
trace-square candidates +-phi and +-(phi-1) violate one embedding.)
Then `v^2=d`. Reflection integrality on unimodular M forces `v in dM`,
so its norm lies in `d^2 O`, impossible for nonunit d. Thus every
crossing is orthogonal. Conversely each primitive M reflector has
norm dividing2 and extends to R. This proves the full arrangement
identity at the actual coordinate facet by local finiteness.

In particular, under the SAME literal M_8 cutoff, every positive
one-coordinate extension `M_8+<d>` with nonunit d is excluded at once.
Changing d, its prime factors or its size cannot evade that barrier.
This statement does not classify nonunimodular positive root modules,
nonorthogonal extensions, or arbitrary lattices in the same K-space.
