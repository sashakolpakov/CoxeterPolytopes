# Hostile L4 audit of the high-rank complete K3-root transfer

**Frozen scope date:** 2026-08-09.  
**Input under audit:** `K3Root[11,20]`, the class of complete all-`(-2)`
root chambers in even rank-11--20 hyperbolic lattices admitting a primitive
K3-lattice embedding.  
**Verdict:** the mathematical 19-row classification is valid, and the marked
K3 carrier/return theorem is valid, but the package is **not strict L4** under
the refined `P -> all carriers -> independent AG -> theorem on P` contract.

## Executive verdict

The high-rank package closes three substantial inverse layers.

1. Every arithmetic input `(N,C)` in the declared class has an exact-Picard
   projective K3 carrier.
2. A Weyl change of the Néron--Severi marking returns exactly the same complete
   chamber `C`; because `C` is cut by every square-minus-two root, there is no
   hidden-root ambiguity.
3. For the nineteen finite-volume rows, primitive embeddings form one orbit
   under both `O(Lambda_K3)` and `O^+(Lambda_K3)`, the period strata are
   nonempty, and independent literal Vinberg audits certify finite volume and
   noncompactness.

The terminal leverage arrow nevertheless fails.  The list of nineteen
lattices and the absence of rank 20 are not derived by an independent
algebro-geometric classification.  Artebani--Hausen--Laface first identify
polyhedral effective cone, finite `O(N)/W_2(N)`, and finite automorphism group;
their next theorem obtains the finite list by importing the already completed
Nikulin--Vinberg classification of the hyperbolic lattices with finite
`O(N)/W_2(N)`.  Roulleau's atlas explicitly starts from the same inherited
118-lattice list and then constructs geometric models and complete curve
configurations.  See the primary [AHL paper](https://arxiv.org/abs/0901.0369)
and [Roulleau atlas](https://arxiv.org/abs/2003.08985).

Thus the AG stage supplies a valuable dictionary and genuine geometry, but it
does not independently enumerate the finite-volume locus.  The conclusion
“nineteen in `H10`--`H18`, none in `H19`” is the reflective-lattice theorem
translated into surface language and then translated back.  That is a strong
calibration, not the new AG pushback demanded by L4.

## Gate-by-gate score

| gate | result | exact reason |
|---|---|---|
| `C` | **pass, narrow** | `K3Root[11,20]` is a non-singleton arithmetic lattice/chamber class and does not mention an actual surface, although primitive embeddability in the K3 lattice is deliberately built into membership. |
| `Car` | **pass for the declared K3 species** | Primitive embedding and period surjectivity produce an exact-Picard K3 carrier for every ambient member.  The nineteen terminal rows have a unique oriented primitive-embedding orbit. |
| `Return` | **pass** | K3 Weyl transitivity sends the nef cone to the prescribed integrally marked complete root chamber; the all-root definition closes the hidden-`(-2)`-wall gate. |
| `AG` | **fail: independent exhaustion** | Period/Torelli, the finite-automorphism dictionary, Cox-ring finite generation, and Roulleau's geometric models are AG.  The actual `19/0` exhaustion is imported from the Nikulin--Vinberg `2`-reflective lattice classification. |
| `Pushback` | **fail: new theorem** | The no-`H19` and nineteen-type conclusions are already contained in that input classification.  Literal Vinberg audits verify the listed positives but cannot prove the list exhaustive. |

The weakest gate controls the promotion: this is not L4.  This verdict does
not dispute the upstream classification theorem or any of the nineteen
literal volume certificates.

## The finite-volume iff must be stated correctly

Raw K3-carrier nonemptiness is automatic throughout `K3Root[11,20]`, so it
cannot characterize finite volume.  The correct statement is

```text
C has finite hyperbolic volume
  iff some/every exact-Picard full-nef K3 carrier has finite Aut
  iff Eff(X) is rational polyhedral
  iff Cox(X) is finitely generated.
```

There is an explicit counterexample to the unqualified shorthand.  Let

\[
N_0=U\oplus E_7(-1)\oplus A_2(-1).
\]

It is even of signature `(1,10)` and determinant `6`.  It embeds primitively
in `3U+2E8(-1)`: use one `U`, the coordinate `E7` subdiagram on simple-root
nodes `[0,1,2,3,4,5,7]` of the first `E8`, and the coordinate `A2` on nodes
`[0,1]` of the second.  A generic period therefore gives an exact-Picard K3,
and Weyl re-marking realizes any complete root chamber of `N_0` as its nef
cone.  But determinant `6` excludes `N_0` from all four rank-11 finite-volume
rows, whose determinants are `2,8,32,128`; the accepted Nikulin--Vinberg
exhaustion therefore makes its complete chamber infinite-volume.

So

```text
unqualified K3 carrier exists       = true
finite-Aut/polyhedral/Cox-finite carrier exists = false
finite-volume chamber               = false.
```

This counterexample does not weaken the decorated finite-volume dictionary;
it identifies the necessary decoration.

## Essential surjectivity, orientation, and exact return

Within the declared category of full-nef K3 carriers, essential surjectivity
is genuine rather than per-example recognition.  Starting from any `(N,C)`:

1. the primitive embedding in `Lambda_K3` is part of the arithmetic input;
2. a period away from the Noether--Lefschetz hyperplanes has exact Picard
   lattice `N`;
3. the K3 ample-cone theorem identifies a nef chamber of the complete
   square-minus-two arrangement; and
4. a Weyl element changes the marking so that this chamber is exactly `C`.

For the nineteen finite-volume rows, the stable discriminant-length
inequality and the complement-`U` orientation correction give exactly one
primitive embedding orbit under `O^+(Lambda_K3)`.  The exact-Picard period
dimensions are `20-rho`, from 9 down to 1.  These facts close the marked
carrier and return questions.  They do not make the later enumeration
AG-native.

The return is an equality of integrally marked lattice/chamber data.  It is
not a bijection between surfaces and abstract polytopes: an entire period
locus returns the same marked chamber, and forgetting the marking is
many-to-one.

## Compact and noncompact statements are separate

The compact finite-volume part of this high-rank class is empty before the
nineteen-row list is consulted.  Every `N` has rational dimension at least
eleven and indefinite signature, so Meyer gives a rational isotropic vector.
For a finite-volume complete K3 root chamber, Weyl reduction produces a nef
primitive isotropic class and hence a genus-one-fibration cusp.  Therefore

\[
K3RootFV^{\rm compact}[11,20]=\varnothing.
\]

This is a clean class-wide theorem, but its decisive input is arithmetic
Meyer isotropy rather than independent AG classification, and it is not new.

The finite-volume noncompact locus is the published nineteen-type list in
hyperbolic dimensions 10 through 18.  All nineteen lattices visibly contain
`U`; the literal replays independently contain 1,156 ideal vertices in total
and no hyperideal rays.  The absence of a noncompact `H19` row, unlike
high-rank compact emptiness, uses the inherited reflective-lattice
classification.

## Why the source dependency is terminal for this promotion

The source chain is exact:

```text
Nikulin--Vinberg classify O(N)/W_2(N)-finite hyperbolic lattices
  -> AHL identify that lattice condition with polyhedral Eff / finite Aut
  -> Roulleau geometrically realizes and studies the inherited surface list
  -> the project returns the same complete root chambers and audits volume.
```

The middle arrows are meaningful algebraic geometry.  But the only arrow
that supplies the finite list and the rank-20 zero is the first, already
polytope/lattice-theoretic arrow.  Consequently the round trip has not
extracted a new restriction from the algebra of the surfaces.

Neither of the following repairs this:

- verifying every positive row by Vinberg's endpoint and parabolic-
  neighbourhood tests; or
- writing the reflective condition as finite generation of the Cox ring.

The first certifies positive membership but not exhaustion.  The second is a
correct equivalence, not an independent classification proof.

## Exact missing gate

An honest L4 promotion needs one of the following.

1. An independent geometric classification of the relevant exact-Picard K3
   loci—using surface models, genus-one fibrations and Mordell--Weil groups,
   curve configurations, Cox algebras, or another genuinely AG mechanism—
   which derives the finite list without invoking the Nikulin--Vinberg
   `2`-reflective lattice exhaustion.
2. A broader raw polytope class for which the K3 transition is still
   exhaustive, but AG effectivity or moduli conditions remove or construct
   cases not already fixed by the Coxeter/lattice input.
3. A new dimension/facet restriction or a genuinely new terminal positive
   realization that appears only after this independent AG stage.

If an independent geometric proof merely reproduces the known nineteen rows,
it would be a valuable alternative theorem, but the strongest version of the
user's payoff gate still asks for a new conclusion about polytopes.  The most
promising promotion therefore enlarges the carrier alphabet beyond complete
`(-2)` K3 nef chambers—mixed walls, invariant sections, or another variety
type—while retaining the exact return theorem.

## Novelty boundary

The honest project-level achievements remain substantial:

- a uniform marked `P -> K3 -> same P` theorem for a non-singleton class;
- a clean AG dictionary for volume, automorphisms, effective cones, and Cox
  rings;
- exact oriented embedding/period bookkeeping for all nineteen rows; and
- literal Vinberg certification of all nineteen positive chambers.

The nineteen abstract polytopes and the no-rank-20 classification are
published known results.  This audit finds no new abstract polytope and no new
AG-derived dimension bound.

## Replayable evidence

- [hostile-audit manifest](../data/enumerations/k3_high_rank_l4_promotion_hostile_audit.json)
- [certifier](../scripts/certify_k3_high_rank_l4_promotion_hostile_audit.py)
- [tests](../tests/test_k3_high_rank_l4_promotion_hostile_audit.py)
- [upstream classification theorem](k3_high_rank_full_root_classification_transfer.md)
- [upstream oriented arithmetic theorem](k3_high_rank_oriented_monodromy.md)

Run:

```bash
python3 scripts/certify_k3_high_rank_l4_promotion_hostile_audit.py --check-archive
PYTHONPATH=src:. python3 -m pytest -q tests/test_k3_high_rank_l4_promotion_hostile_audit.py
```

The optional `--check-sources` flag rechecks the pinned local AHL and
Roulleau source copies and exact dependency windows.
