# Root review: mixed right-angled carrier errors and maximal ideal degree

19 October 2026. **PASS, full written-source review.** No mathematical
enumeration, optimizer, Lean module or typesetting was used for this review.
The two sources reviewed in full are:

* `right_angled_mixed_nu12_carrier_error_sharpening_research_20261019.md`,
  SHA256 `e638488168447f03f9619935aede2e6854523b147e33c855b785e9f52e83eed8`;
* `right_angled_mixed_nu12_carrier_forced_by_maximal_ideal_degree_20261019.md`,
  SHA256 `7c6d300a5f2df7ba0d1f1f7ae1154b9b45c7e8163c83347c8efeb8c90413b233`.

The earlier local source `fc047b78e9a8fc9eaf283bd5d689b69d72a8bbf0976a4765600ffbe8625955b0`
already has a separate root review. I reread its complete argument and
the actual P3/quad cross identities, and checked the quoted triangle
density constants against the existing weighted-face proof.

## Universal error bounds

The identities are `p+4c=2nu-4`, `4c=E_FI+2E_II` and
`3p=E_FI+2E_FF`. For induced-path finite corners, `E_FI>=3Z`,
`E_II>=2Z`, and `p>=Z` give `4Z<=nu-2`. The attachment lemma is
essential: a bad ideal/ideal edge cannot be reused outside P3.

At nu12, two path corners force c4,p4. Their four distinct bad
ideal/ideal edges give E_FI<=8, while their own six incident edges
and the remaining two finite vertices give E_FI>=10. Thus Z<=1.
The integral nu13 case and the real nu>=14 bound give
`Z<=nu-11` throughout the remaining range.

For the companion count A, `3A<=4c`, `p>=A` give
`A<=(nu-2)/2`. At nu12 or13 the only exceptional cusp counts are
three or four. Three cusps and three such vertices give an actual
K3,3 subgraph of the planar edge graph. At nu12 with four cusps,
three all-ideal-edge finite vertices force all four to have that
property: there is no available finite neighbor for the fourth.
There are then only two ideal/ideal edges but at least four triangular
facets; all would be bad, contradicting the attachment lemma. The
nu13 four-cusp case similarly forces six such vertices and eighteen
finite/ideal edges, exceeding capacity sixteen. This verifies the
small endpoints and hence `A<=2(nu-11)`.

The summations consequently give `sum M<=E12` and companion count
at most `2E12`, without finite-edge errors. They change the error
allowance, not the definition or multiplicity of an actual quad move.

## Maximal ideal degree and the forced carrier

Under the temporary assumption that every three-face has nu<=11,
the local theorem makes the graph on the ideal directions complete
minus a matching. A P3 move preserves all ideal directions: an
untouched one has at least two preserved cross bad edges to its
triple. A missing-pair quad preserves all the ideal directions by
its exact cross identities. Taking a globally maximal ideal degree K
before making any move rules out additional ideal directions at the
opposite vertex. Thus the class of K-maximal vertices is closed.
This maximality step is necessary; no constancy of finite directions
is assumed without it.

I checked the displacement coefficients directly by removing triples
through matching edges and adding one half-weighted move per actual
missing pair. The spatial coefficients outside the K ideal directions
are zero. For K9--12 the stated four causal margins and their integer
minima agree with this calculation. The K9,r0 endpoint is a nonzero
future null vector and is sufficient for the cancellation argument.
Undirected weights agree at opposite endpoints because the ideal set
and the matching-complement property persist there. A finite closed
component cannot have a zero sum of nonzero future causal vectors.
It follows that K<=8.

Finally `2B=sum t(G)<=binom(8,3)f=56f` contradicts the previously
reviewed lower bounds `7545f/118` in dimension11 and `5030f/59`
in dimension12. The finite-vertex count is positive in these
dimensions, since an entirely ideal example is excluded already by
the published dimension-six ceiling. Therefore every hypothetical
dimension11/12 example contains a three-face with nu>=12.

## Verification boundary

These are reviewed written necessary-condition theorems. They do not
exclude dimension11 or12, bound the density of the required carriers,
or discharge the remaining unrestricted incidence/Gram problem.
They have not been integrated into the primary mixed manuscript or
formalized in Lean. The unconditional ceiling remains12.
