# Independent audit of the D-center common-opposite estimate

20 October 2026. Reviewer: resource_cleanup. **PASS**, within the stated
compact, actual-wall scope. This is a full written mathematical audit,
not an execution, numerical test, or formal verification. No source or
global status was changed.

## Exact source binding and inputs

Reviewed in full:

* `results/compact_dn_center_neighbor_common_opposite_energy_20261020.md`,
  SHA256 `e371611c4903ad362654ee0d4656457140cdc5e4374b4e19a2697718c2d679b1`;
* its energy/pair dependency
  `results/compact_universal_e6_e7_e8_component_energy_20261020.md`,
  SHA256 `acf261d97f6eb1a11a64502e17e079ef7f418b9746fbf85dbdf6915750fe9a18`.

I also read the underlying actual reciprocal-shadow and farthest-energy
proofs. The classical finite Coxeter classification, exclusion of actual
affine principal diagrams in a compact polytope, and Lanner order at most
five remain explicit external inputs. This audit does not replace their
classification proofs. No arithmetic or field restriction is used here.

## Complete classification, including spectators

The positive definite replacements deleting `a` and deleting `o` are
essential. Every old node belongs to at least one of them, so every
contact considered is finite. The first replacement leaves the right
`A3`; its only possible nonempty contact is middle3 or end3/4/5. The
second permits at most one left-chain contact. Cycles and the
4/3/3 star give nonpositive principal blocks in these finite
replacements, rather than merely indefinite diagrams in the full star.

I checked each remaining branch independently:

* Middle3 forces contact4 at `a`; an extra spectator produces the retained
  affine-F4 minimum path `3,4,3,3`.
* End5 gives a retained H4, hence no other white contact. A contact at
  `a` produces the retained `3,3,3,5` path, whose raw determinant is
  `6-4phi<0`.
* End4 permits no additional retained contact, forces a contact at `a`
  by the exact affine terminal-4 fork, and that contact must be3.
* End3 excludes all left contacts. A spectator3 gives the full affine
  star with arms1,3,3; a spectator4 gives the terminal-4 fork. Without
  either, isolated3 has negative height square at rank7, and contains
  affine E8 from rank8 onward.
* With no right contact, all left contacts except an end3 are excluded
  by retained finite obstructions or exact affine forks. In the last
  case a minimal nonpositive subdiagram must contain both `c,o`, since
  deletion of `o` is positive definite. Its connecting path already
  has at least six vertices. Spectators attach through `c` and cannot
  shorten it, contradicting the imported Lanner bound.

Thus A/B/C exhaust the complete old contacts. The monotonicity argument
is valid where invoked: the connected base obstruction has a positive
nonpositive quadratic witness. Extra contacts do not repair it.

## Edge pressure and the D7 conclusion

Using the unit inverse-D entries, I obtain exactly
`delta_A=1/k`, `delta_C=1/(2k)` and, in both cases,
`cosh L_o=a+k/sqrt(2)>7`. Hence their **half**-tangent squares exceed
`3/4`. For isolated5 I obtain

    tanh^2 L_a = 1-(9-4sqrt(5))/(r-3).

At rank7 its half-tangent square is `4/sqrt(5)-1>3/4`; the stated
`256>245` comparison is correct. Since covering energy is at least
each individual half-tangent square, every common-opposite profile
has energy greater than `3/4`.

For D7 the dual correlation is `k=sqrt(4/5)`. At energy at most `1/2`,
the sharper pair inequality gives
`tanh L_o*tanh L_a <=16(1+k)/(5+k)^2<k`. Clearing denominators gives
`k^3+10k^2+9k-16=(49/5)k-8>0`, with the stated exact square seal.
The actual reciprocal shadow therefore forces a common opposite,
contradicting the preceding classification pressure. **E_D7>1/2 is
proved with arbitrary finite spectators.**

The D8 boundary is correctly left open. At its rank/14 target `e=4/7`
and `k=sqrt(5/6)`, the same pair upper bound is
`56(1+k)/(9+2k)^2`; the cleared difference from `k` is
`(85/3)k-26<0`, since `36125<36504`. Consequently this estimate does
not force the two opposites to coincide. Neither D8 at that target nor
an unrestricted compact dimension ceiling follows from this source.
