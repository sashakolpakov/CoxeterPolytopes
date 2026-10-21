# The long-end isolated-eight star at an arithmetic B vertex

20 October 2026. Own written composition, **awaiting outside review**.
No computation, enumeration, build or formalization. This excludes a
specified vertex star, not a dimension range of arbitrary polytopes.

**Claim.** A compact arithmetic Coxeter n-polytope of first type, n>=9,
cannot have a whole finite B_n vertex and an actual outside facet whose
only nonzero contact with that vertex diagram is a finite label8 at
the long-chain endpoint. In the standard B_n simple roots
`beta_i=-e_i+e_(i+1)` (i<n), `beta_n=-e_n`, this is node1.
Other isolated-eight contact positions are not included.

The two substantial inputs are the complete invariant-lattice proof

* `silver_b_vertex_invariant_lattice_flags_and_even_neighbor_20261020.md`,
  SHA256 `2fed0220943f80acd6747119e98e1edee39d2ea4fe9073acff5ce6b894770264`,

and the reviewed literal silver-family obstruction

* `odd_silver_stable_integral_isometry_and_uniform_obstruction_20261020.md`,
  SHA256 `2ad82b9a52c7b81364a070b985932932eab72438922e7b0d717b392e5e4c5ee3`.

All paths are under `results/`. The latter includes the explicit
Hasse/strong-approximation inputs and the earlier closed E8-host
geometry. They are written mathematical inputs, not executions made
for this composition.

## 1. The actual star fixes the field and quadratic space

Write p=sqrt(2), alpha=1+p, X=p*alpha=2+p. The cyclic product at
the label8 bond is `4cos^2(pi/8)=X`, so the arithmetic ground field
k contains K=Q(p). Work in its admissible k-quadratic model.

The old long roots can be scaled over k to have norm2 along their
3-chain. A totally positive scalar change of the form performs the
initial common normalization. At the terminal4-bond, a short root of
norm b and contact g satisfies `g^2=b` after the long root has norm2;
multiplication of the short root by `-1/g` therefore gives norm1 and
contact -1. Thus the
whole old B_n is literally the signed-permutation root system on an
orthogonal unit basis e1,...,en defined over k.

Let an initial cap root have norm q and contact g with beta1. The
native8 equation is `2g^2/q=X`. Rescale it over k to make its contact
`-alpha`; its norm becomes `2alpha^2/X=X`. Orthogonality to the
other old roots makes its spatial component `alpha e1`. Subtract
that component and call the remaining vector e0. It is orthogonal
to every ei and has norm `X-alpha^2=-alpha`. The entire full-rank
star has therefore produced the explicit k-basis and form

    <-alpha> + I_n,
    c=e0+alpha e1.                                    (1)

Its Gram is over K and has Lorentz signature at the chosen place.
If k properly extended K, a nonidentity real embedding of k fixing
this copy of K would leave the full-rank Gram unchanged. This
contradicts positive definiteness of an admissible form at every
nonphysical embedding. Consequently k=K. This is the same full-rank
cyclic-field argument used in the separately reviewed isolated A_n
packing proof; no arithmeticity is inferred for an arbitrary group.

Spatial sign changes turn the cap in (1) into `e0-alpha e1`.
Conjugating by all old signed permutations gives exactly the W0
of the invariant-lattice input. This is a full-rank actual-wall
normalization, not just a match of one edge length or determinant.

## 2. The arithmetic group preserves an actual lattice

Let Gamma be the polytope's arithmetic reflection group. In the basis
(1) it is commensurable with the integral orthogonal group of some
full O_K lattice L0. There is an invariant full lattice for Gamma:
take H=Gamma intersect O(L0), which has finite index in Gamma, choose
left-coset representatives gamma_i, and put

    Lambda=sum_i gamma_i L0.

This is a finitely generated full O_K module and is Gamma invariant,
since Gamma permutes the summands up to H. In particular it is W0
invariant. No automatic reflectivity transfer between commensurable
lattices has been assumed.

The flag classification now gives, in the SAME coordinates,

    Gamma subset O(L), L=<-alpha>+I_n,

or, only if n is even,

    Gamma subset O(S), S=L+O(0;1,...,1)/p.             (2)

Multiplying Lambda by a scalar in that proof does not change its
orthogonal group. Every two-level dyadic lift is covered; assuming
that the original lattice was L would leave a real gap here.

## 3. Compactness contradicts the two full arrangements

All actual reflection mirrors of Gamma belong to the full integral
reflection arrangement of the group on the relevant right side of
(2). A full-arrangement chamber meeting the interior of the compact
Gamma chamber is contained in that compact chamber. The integral
group is discrete and its mirrors are locally finite, so this full
chamber is compact as well.

For O(L), the inherited silver-family theorem gives infinite-volume
full chambers for every n>=9, a contradiction. For O(S), n is even,
so n>=10. The input's even-neighbor proof chooses an actual coordinate
facet and proves that **every induced mirror is a standard L_(n-1)
mirror**. Hence the facet contains an entire standard chamber, whose
volume is infinite for n-1>=9. Such a facet cannot be compact. This
again contradicts compactness of the full chamber.

This proves the stated claim subject to the full outside review of
the flag theorem and composition. It is stronger than the earlier
literal-lattice obstruction: it handles every invariant lattice
compatible with this actual star. It does not exclude isolated-eight
caps at nodes2 or3, other cap labels or types, all silver forms, or
an entire rational commensurability class. It changes neither the
compact29 nor the arithmetic compact20 ceiling. No new manuscript
or Lean package is proposed for this local result.
