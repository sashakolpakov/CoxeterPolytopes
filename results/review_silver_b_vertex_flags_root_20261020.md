# Root review: silver B-vertex flags and the even neighbor

20 October 2026. **PASS, full written review.** No enumeration, script,
build or Lean execution. This is an invariant-lattice classification
and a specified-star obstruction; it gives no uniform dimension bound.

Reviewed source:
`results/silver_b_vertex_invariant_lattice_flags_and_even_neighbor_20261020.md`,
SHA256 `2fed0220943f80acd6747119e98e1edee39d2ea4fe9073acff5ce6b894770264`.
I read the complete candidate before the congruence clarification and
checked the corrected formula afterward. The bars now correctly mean
top residues modulo p with canonical 0/1 lifts, not residues modulo2.
The recovery criterion and subsequent proof were unchanged.

## Lattice sandwich, all flags and lifts

At odd primes the sign changes split the coordinate ideals, permutations
identify the spatial ideals, and the cap's two unit off-diagonal entries
identify the time ideal. PID normalization removes these local scalars.
At the dyadic prime, `p^2E00` follows from `p^2I-sum p^2Eii`.
Multiplying by cap differences gives all `p^3` matrix units. A primitive
vector first yields `p^3L`, then the stated spatial-unit and time-only
cases give every `p^2ei`. The two-level sandwich is proved without an
unbounded module search.

Because `O/(p^2)=F2[p]/(p^2)`, the displayed top/middle/lift description
is exact. Any nonconstant spatial code vector generates the whole
augmentation module by transpositions; otherwise it lies in the
two-dimensional fixed space. These give exactly the ten listed codes.
If the top code is not D, the cap images and permutations force H into
the middle code. I checked every allowable middle code, every top code
inside H and every equivariant lift. In particular the pair-vector
relation kills `Hom(E,F2)`, and transposition differences cannot be D
or the all-ones vector when n>=3. No lift family is omitted.

## Canonical recovery, rather than an isometry assumption

For integral y, expansion modulo p^2 gives
`Q(y)=sum bar(x_i)-(1+p)bar(t)` modulo p^2. Hence divisibility by
p^2 is precisely the code E. This proves the first-level recovery
formula. I checked all seven rows of the dual/norm-hull table,
including the parity-dependent recovery of S from J and Z.

For top D the fractional generator has norm valuation -2 or -1;
its cross term with an integral vector is integral. The integral-norm
set of the scaled lattice is therefore exactly L_T, for every lift.
For top H, the dual's deepest generator has norm valuation -4 or -3,
and its once-divisible version still has negative valuation. The
integral-norm hull is L. For top E, that deepest term cannot cancel
against the other norm or cross terms, so the remaining first-level
Z hull gives L or S. In the final even-J case I rederived the dual
condition `sum y_i=a*b mod p`; integral norms force b divisible by p
and give exactly K0. Its dual is literally L_Z, recovering S.

These operations are invariant under the FULL orthogonal group and
are genuine lattices in every use. At odd primes they recover L,
so they prove global group containment in the stated coordinates.
They do not merely identify a rational quadratic space.

## Even-neighbor facets

The primitive S-reflector pool 1,2,2+sqrt(2) follows from
`2v/q in S# subset L subset S`. A fractional norm2 reflector is
impossible. Bounding both embeddings gives exactly the integer and
half-integer coordinate lists; I checked every projection scaling.
The result is a SUBSET of standard diagonal mirrors, which is the
needed direction: an induced facet contains a full standard chamber.
No extension of every intrinsic reflection is needed here.

The zero-height trace lower bound n excludes all fractional roots
for even n>=10. Integral norm2 coordinate multiples sqrt(2)*ei
are imprimitive and add no mirrors; the full zero arrangement is B_n.
Thus a coordinate mirror is an actual facet. If the ambient chamber
were compact, its facet would contain an infinite-volume standard
L_(n-1) chamber by the explicit inherited theorem, a contradiction.
No intrinsic Coxeter property of that metric facet is assumed.

The proposed long-end composition is also correct in mathematical
scope: the full-rank actual native8 star fixes K and the displayed
quadratic space; arithmetic commensurability supplies an invariant
lattice by a finite LEFT-coset sum; the flag theorem applies to that
lattice. Inclusion of actual reflection mirrors makes a full integral
chamber lie inside the compact subgroup chamber. This supplies the
contradiction without treating reflectivity as commensurability
invariant. Its source and separate outside review bind the final
composition bytes. Nodes2 and3 and all other labels require their
own arguments; none is included in this review's geometric scope.
