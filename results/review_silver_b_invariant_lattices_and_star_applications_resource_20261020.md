# Outside-author review: silver B lattices and the node-one/node-three stars

**PASS, subject to the explicitly inherited literal silver-family theorem.**
I fully read and independently rederived the three new sources below.
I found no mathematical gap after the stated top-residue clarification.
This is written mathematical review, not a scientific execution or a
new uniform compact dimension bound.

|Source under `results/`|SHA256|
|---|---|
|`silver_b_vertex_invariant_lattice_flags_and_even_neighbor_20261020.md`|`2fed0220943f80acd6747119e98e1edee39d2ea4fe9073acff5ce6b894770264`|
|`silver_long_end_eight_arithmetic_vertex_exclusion_20261020.md`|`cc7f93302e7a81eda7e7b2dd8fbc88b79b923958fab7bde014cbe43a6bd34546`|
|`silver_b_node_three_triple_cap_lattice_extension_20261020.md`|`47ec9726572e9ccc5e87295689849e45a2b23c2499a3457beaa4cd2bffb6ce8f`|

I also read the inherited scope and proof of
`odd_silver_stable_integral_isometry_and_uniform_obstruction_20261020.md`,
SHA256 `2ad82b9a52c7b81364a070b985932932eab72438922e7b0d717b392e5e4c5ee3`.
Its genus/strong-approximation and closed E8-host inputs remain inherited
written theorems; this review is not a fresh independent verification
of their entire recursive proof chain or primary-source provenance.

## 1. All invariant lattices, including the two-level lifts

At an odd prime, the signed spatial permutations split the coordinate
ideals and identify all spatial ones. The singleton reflection has unit
off-diagonal entries there, identifying time with space. A global scalar
can remove these odd-prime ideals over the PID. This changes no
orthogonal group.

At the dyadic prime, `p^2 E00` follows from `p^2 I` and the spatial
matrix units. Multiplication by the cap difference yields the claimed
`p^3 E0i,p^3 Ei0`; hence primitivity supplies `p^3 L`. I checked both
primitive-coordinate cases improving this to `p^2 L`: the unit spatial
coordinate case and the time-only unit case require different uses of
the cap, and both are present. No case has been omitted.

Over `O/(p^2)=F2[p]/(p^2)`, the description by `U subset T` and an
equivariant lift `f:U -> V/T` is complete. Sign changes are identity
modulo p^2. The cap condition is precisely `Ni U subset T`, independent
of the lift. The ten invariant codes follow by taking transposition
differences and then passing either to `V/E` or to the constant-spatial
subspace Z. For U other than D, a nonzero cap image and permutation
invariance force `T superset H`.

I checked every displayed lift: the two H lifts; the zero lift on E;
the two even-n J lifts; and all five possible middle codes above D.
In particular the transposition argument for fixed classes modulo D/Z
needs n at least three, exactly as stated. The pair-vector relation
correctly gives `Hom(E,F2)=0` even in characteristic two.

The final source now defines the bars in the norm congruence as residues
modulo p, canonically lifted to 0 or 1. This corrects the earlier literal
modulo-2 wording. With that definition,

`Q(y) in p^2 O iff y mod p belongs to E`

is exact, and all uses of the integral-norm hull follow.

I checked every row of the first-level recovery table, including the
bilinear duals. For the remaining lifts:

* U=D: `Q(d)=n-alpha` has valuation zero for even n and one for odd n.
  The fractional generator in `p^-1 Lambda` therefore has negative norm
  valuation for a unit coefficient; its integral cross terms cannot
  cancel it. The hull is exactly L_T for every lift.
* U=H,T=H: the two congruence kernels and their dual generators are
  correct. Coefficient valuations zero or one both give negative norm
  valuation, leaving exactly L.
* U=E,T=H: a unit coefficient of `w/p^2` has valuation -4 or -3,
  whereas the other square and cross terms have valuations at least
  -2 and -1. The surviving first-level module is `p^-1 L_Z`.
* U=J,T=H, n even: direct pairing gives exactly equation (9).
  After its leading coefficient is divisible by p, the norm criterion
  gives `K0=O e0 + p^-1{spatial sum in pO}` for both lifts. Its dual
  is L_Z, whose canonical recovery is S.

All these hulls are explicitly lattices, so the use of the integral-norm
set as a canonical operation is legitimate. The odd-prime normalization
and the dyadic recovery combine by intersection of local lattices. This
proves literal same-coordinate containment of the full isometry groups
in O(L), or O(S) only when n is even. It is neither homothety of all
lattices nor a reflectivity inference from commensurability.

## 2. The even-neighbor facet argument

For even n, Q(S) is integral, although its bilinear pairings need not
all be integral. The source correctly uses the bilinear dual:
`S# subset L subset S`. Primitive reflection integrality therefore
forces q to divide 2, giving the normalized positive pool `1,p alpha,2`.
A fractional norm-two reflector would have to lie in S# and hence L,
which excludes that case.

Both real strict inequalities for a mirror crossing the coordinate
hyperplane are justified by independence and positive definiteness at
the other embedding. Summing them gives the stated trace bounds. I
checked the complete integral and fractional k tables directly; the
conjugate inequality is necessary to remove the opposite-sign mixed
half-coordinate cases. Every projected ray has the stated integral
scaling and resulting norm. In particular `2/(p alpha)=2-p` is integral
in the unscaled integral case. Thus each induced mirror is indeed a
standard L_(n-1) mirror.

At the time axis, fractional zero-height norms have trace at least n.
For n at least ten this excludes them all. Integral zero-height roots
give only the B_n mirrors; a possible imprimitive vector p e_i adds
no new mirror. Hence the chosen coordinate mirror is an actual facet
of a full-arrangement chamber.

The inclusion direction in the final argument is correct and important:
the facet's bounding mirrors form a SUBSET of the standard arrangement.
A generic point of the facet therefore determines a whole standard
chamber contained in that facet. Noncrossing ambient mirrors have
constant sign on the coordinate hyperplane and add no hidden boundary.
The inherited infinite-volume standard chamber cannot fit in a compact
facet. Neither surjectivity of mirrors nor an intrinsic Coxeter-facet
assertion is required.

## 3. The actual long-end star and its arithmetic lattice

The old 3-chain is normalized by root scalings over the ground field
and a totally positive common scalar change of the form. The terminal
short root is multiplied by `-1/g`, yielding norm one and contact -1.
This gives a literal orthonormal spatial basis over that field.

For node one, the native-eight equation is `2g^2/q=X=2+p`.
The prescribed scaling gives cap norm X, spatial component alpha e1
and residual time norm `X-alpha^2=-alpha`. The star spans the full
quadratic space. A proper totally real extension of Q(p) would have a
nonphysical embedding fixing this full Lorentz Gram, contradicting
admissibility. Thus the asserted ground field and basis are justified.

The finite sum of left-coset translates `sum gamma_i L0` is a genuinely
invariant lattice: left multiplication permutes the representatives
up to the stabilizing subgroup on the right. No normality assumption
on that subgroup is needed. The flag theorem therefore applies to the
actual arithmetic group rather than to a guessed integral model.

Containment of reflection groups in O(L) or O(S) puts a full integral
arrangement chamber inside the compact actual chamber. Local finiteness
makes its closure compact. The inherited L-family obstruction, or the
even-neighbor facet argument, gives precisely the claimed contradiction.

## 4. Node-three transfer

For n at least four, triple indicators span the spatial module over the
dyadic DVR: exchanges give all augmentation differences and the sum
three is a unit. This also proves the residue-field span used in the
primitive spatial-coordinate case. The matrix-unit construction and
both primitive sandwich cases are valid.

Modulo p^2 the triple cap has image `(t+sum_T xi)(e0+j_T)`. Its common
kernel is D, and its image orbit spans H. Consequently it has exactly
the same flag criterion as the singleton cap, independent of the lift.
Every triple-cap invariant lattice is therefore also singleton-cap
invariant at the dyadic prime; the odd-prime scalar lattice is already
preserved. This is an exact global lattice-preservation conclusion.

The node-three normalization yields cap norm `2-p`, spatial vector
`e1+e2+e3`, and again time norm `-alpha`. The full-rank field argument,
invariant-lattice construction and compactness contradiction then apply
as written. Node two instead has time norm `-p` and is not covered.

**Limits:** the reviewed applications concern only isolated native-eight
contacts at nodes one and three of a whole B_n vertex, n at least nine,
under arithmetic first-type hypotheses. They establish neither a uniform
compact dimension ceiling nor a theorem about all silver forms, labels
or commensurability classes. No candidate source was changed by me and
no executable import, scientific run, build or formalization was made.
