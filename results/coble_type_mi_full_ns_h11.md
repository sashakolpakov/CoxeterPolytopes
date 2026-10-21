# Type MI gives the Fermat `H11/206` chamber as a complete Coble nef cone

**Status (2026-08-09): complete positive surface-to-polytope theorem, exact
marked identification with an existing project polytope, and independent GMP
double description.**  This closes the Type-MI row of the nine-class Coble
program.  It is a new complete realization in the project, not a new abstract
polytope or arithmetic group; worldwide priority for this Coble interpretation
remains open.

## Theorem

Let `S` be Kondō's unique Type-MI strong Coble surface in characteristic
three.  Then the projectivized full Néron--Severi nef cone of `S` is a
finite-volume noncompact arithmetic Coxeter polytope

```text
P_MI in H11,       facets(P_MI)=206.
```

Its complete negative-curve split is

```text
2*(-4) + 30*(-2) + 174*(-1).
```

Moreover an explicit facet permutation and a determinant-one integral
isometry identify `P_MI`, with orientation and every wall retained, with the
independently constructed Fermat-quartic chamber `F3-C2-H11-206`.  Thus the
two constructions give the same Coxeter polytope, not merely lattices with
the same determinant or diagrams with the same facet count.

The complete [certificate](../data/enumerations/coble_type_mi_full_ns_h11.json)
records all divisor vectors, all 247 Coble--Mukai cusps, the 206-facet map,
and the marked lattice isometry.

## 1. Surface and divisor lattice

Kondō constructs `S` by blowing up the ten intersection points of curves of
bidegrees `(1,3)` and `(3,1)` on `P1 x P1`.  In the basis

```text
f1, f2, e1,...,e10,
f1.f2=1, ei^2=-1,
```

the boundary and canonical classes are

```text
B  = f1+3f2-sum(ei),       B^2=-4,
B' = 3f1+f2-sum(ei),       B'^2=-4,
B.B'=0,                    B+B'=-2K_S.
```

Hence `NS(S)=U+<-1>^10=I_(1,11)` exactly.  The thirty published `(1,1)`
curves are reconstructed without floating-point finite-field calculations:
the fifteen duads and fifteen synthemes on six letters determine thirty
four-subsets of the ten complementary triads, and therefore thirty classes

```text
f1+f2-e_i-e_j-e_k-e_l
```

of square `-2`.  Their complete incidence Gram is Kondō's
duad--syntheme part of `Gamma_MI`.  The ten ordinary exceptional curves meet
both boundaries once.

## 2. The missing walls come from an exact cusp lift

Form the saturated extended Picard lattice

```text
Pic~ = <NS(S), B/2, B'/2>
```

and compute

```text
CM(S)=B^perp intersect B'^perp inside Pic~
```

by a Smith kernel, independently of the forty displayed roots.  The result
is `E10`.  The thirty actual `(-2)` curves together with the ten folded roots

```text
B/2+B'/2+2e_i
```

have rank ten and index two in `CM(S)`.  Exact GMP double description of
their complete chamber gives

```text
CM finite vertices       6,210
CM primitive cusps         247
CM hyperideal rays            0.
```

For each primitive cusp `f` the two possible physical lifts are

```text
E_B  = f-B/2,       E_B' = f-B'/2.
```

The integral census is

```text
165 cusps: neither lift is integral;
 82 cusps: both lifts are integral.
```

Thus the 82 split cusps produce 164 additional `(-1)` curves.  Together with
the ten ordinary exceptional curves, the full list has 174 `(-1)` curves.

This list is exhaustive, rather than a collection of plausible walls.
Kondō proves that the forty displayed classes are all effective irreducible
Coble--Mukai roots.  Every other negative curve must be a `(-1)` curve
meeting one boundary twice; it maps to a primitive nef isotropic class
`f=E+B_i/2`.  Conversely the primitive-isotropic/fibration correspondence
and the Halphen blow-up theorem recover the exceptional component from every
integral `f-B_i/2`.  Hence there are no hidden negative curves.

## 3. Exact polytope and literal Vinberg gates

The 206 physical walls have rank 12 and Lorentzian signature `(1,11)`; their
rectangular wall-system Gram has nullity 194.  They are acute and every intersecting
pair has a Coxeter angle.  The class

```text
h=3f1+3f2-sum(ei),       h^2=8
```

pairs positively with every wall.  Orthogonally projecting `h` to each wall
gives a timelike relative-interior point of that facet, with all other wall
pairings positive.  This is an exact irredundancy witness for all 206
inequalities.

The newly frozen independent GMP ray archive gives

| extreme-ray type | count |
|---|---:|
| finite | 327,566 |
| ideal | 18,354 |
| hyperideal | **0** |
| total | 345,920 |

The [exact-cone certificate](../data/enumerations/ss_fermat_quartic_char3_c2_h11_exact_cone_certificate.json)
and [complete compressed ray archive](../data/enumerations/ss_fermat_quartic_char3_c2_h11_exact_cone_rays.jsonl.gz)
are now frozen.  The latter has SHA-256
`07f7b280aa339bcc4eaadc7b885b1492434b79a01e7464957f4db55546eb06bc`.

The logically independent literal Vinberg pass enumerates all elliptic
rank-ten edge diagrams and all affine zero-neighborhoods:

```text
elliptic rank-10 edges       2,688,957
endpoint types FF/FI/II     1,177,760 / 1,247,706 / 263,491
connected affine components    60,458
bad edge endpoint counts             0
unrecognized ideal endpoints          0
invalid parabolic neighborhoods       0.
```

Every edge therefore has exactly two finite or ideal endpoints, and every
maximal parabolic neighborhood is a full rank-ten cusp neighborhood.  This
is Vinberg's finite-volume criterion in literal form, not an inference from
the absence of an obvious Lannér obstruction.  The ideal vertices make the
polytope noncompact.

## 4. Arithmetic lattice and exact Fermat identity

Boundary reflections do not preserve the odd physical `NS(S)`, as they
should not: a boundary has square `-4` and meets an exceptional curve once.
On the canonical arithmetic lattice

```text
L=2*Pic~ = U + E8(-2) + A1(-2)^2
```

all 206 primitive mirrors reflect integrally.  Exact invariants are

```text
signature       (1,11)
determinant     -1,024
Smith diagonal  1^2 2^10
Nikulin type    (12,10,1).
```

The mirror types are

```text
164 x (square -2, divisibility 2)
 30 x (square -4, divisibility 2)
 12 x (square -2, divisibility 1),
```

and the mirrors span all of `L` by Smith diagonal `1^12`.

These data suggested the Fermat match, but do not prove it.  The certificate
instead orients the independent `F3-C2-H11-206` mirrors, applies an explicit
206-term facet permutation, and checks all `42,436` Gram entries.  From twelve
independent matched mirrors it then constructs a unimodular matrix `T` and
checks

```text
T*(all Type-MI mirrors) = all matched Fermat mirrors,
T^t*G_Fermat*T          = G_(2 Pic~).
```

Kondō's canonical double cover of the Type-MI surface is the supersingular
Fermat quartic, so the exact lattice calculation closes the expected
geometric lineage all the way to the same chamber.

## 5. Novelty and scope

This is a substantial positive result, with deliberately narrow novelty:

- the Type-MI surface, its Fermat cover, and its forty Coble--Mukai roots are
  known;
- `F3-C2-H11-206` was already a certified project polytope;
- its lattice and maximal reflection group are known: Alexeev's
  `(12,10,1)` 14-wall chamber;
- the new result is the complete Type-MI `206`-wall full-nef realization,
  the 247-cusp/164-lift calculation, the exact quotient/K3 same-polytope
  identity, and the missing GMP hyperideal-zero archive;
- no literal 206-wall record occurs in the pinned ordinary or doubling
  CoxIter corpora, but this is only a bounded exclusion.  Worldwide priority
  for the literal subgroup chamber or this surface interpretation remains
  open.

It is still a `SURFACE -> POLYTOPE` theorem (`N/A/L0` on the inverse score).
It is a positive calibration for the proposed `POLYTOPE -> SURFACE` theory,
not the missing non-singleton carrier-universality theorem.

**Directional update.**  The proof in this report remains surface-first.
A separate raw-wall extractor now contracts the intrinsic ten-wall frame,
recovers the Möbius design `S(3,4,10)`, proves the parametric Frobenius
carrier is Coble/root-compatible iff `q=3`, and returns this same 206-wall
chamber.  That companion is a structural P-first `L3` calibration; see the
[remaining-seven audit](coble_remaining_seven_contraction_design_audit.md)
and [independent hostile audit](coble_four_point_constructible_incidence_hostile_audit.md).

## Reproducibility

- [certifier](../scripts/certify_coble_type_mi_full_ns_h11.py)
- [certificate](../data/enumerations/coble_type_mi_full_ns_h11.json)
- [tests](../tests/test_coble_type_mi_full_ns_h11.py)

```bash
PYTHONPATH=src:. python3 scripts/certify_coble_type_mi_full_ns_h11.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_coble_type_mi_full_ns_h11.py
```

Primary source: Shigeyuki Kondō,
[*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591),
especially the Coble--Mukai nef and Halphen results, Example 4.17,
Lemma 4.18, Proposition 4.19, and the Type-MI characteristic/uniqueness
arguments.
