# Hostile audit: PB-Atlas numerical Stage A

**Verdict (2026-08-10): terminal PASS after one repaired API defect.**  The
Schur contraction rule, symmetry quotient, nonsaturated-wall handling,
(\mathbf P^2/\mathbb F_e) endpoint reconstruction, complete coordinate
programs, and all four frozen censuses are correct.  The defect was confined
to an optional expert API: an explicitly supplied list of automorphisms was
checked memberwise but was not required to be a group.  The core now rejects
duplicates and proves generated-subgroup equality before using canonical
masks or stabilizers.  The automatically discovered frozen groups were
already genuine groups, so no calibration count or endpoint changed.

Audited theorem:
[`pb_atlas_numerical_contraction_stage_a.md`](pb_atlas_numerical_contraction_stage_a.md).

## Pinned audited artifacts

| artifact | SHA-256 |
|---|---|
| `results/pb_atlas_numerical_contraction_stage_a.md` | `69bae267db3d9def2887a8aab21acda134a9481e65c3c434ab41be3dcfe148b4` |
| `scripts/pb_atlas_contraction_core.py` | `3787c63a1beabd0b8737da1c66293280b6848fc2761af6abd8ce1a4e57d3107d` |
| `scripts/certify_pb_atlas_numerical_stage_a.py` | `fe4729bb5bc3dcf77a4294893b994c3d85b72c57619577847ddad137eabca2e7` |
| `scripts/freeze_pb_atlas_stage_a_fixtures.py` | `502b823a3f0289a1bef59f0703deda062e50b60267ecfb8d040a26f48a0b870a` |
| `data/fixtures/pb_atlas_numerical_stage_a_inputs.json` | `5b31cb4bf3ef53032f1ce1ff837c2223df40ca6d1cac01a82e42d6ce6aaebd75` |
| `schemas/pb-atlas-physical-branch.schema.json` | `b1bea0891125061c14453f0365993a414afd55aefd9899549f4a5642efc88aca` |
| `schemas/pb-atlas-numerical-stage-a.schema.json` | `6f04b8b116d07874c75125cfb3ed3b791de30cadfa76b90e2c2403c414ec1048` |
| `data/enumerations/pb_atlas_numerical_contraction_stage_a.json` | `6eeb5726cf991e3e7f014429e7a5eb88a7ad4101376b1ced1ef32d30f4b6dcc4` |
| `tests/test_pb_atlas_numerical_contraction_stage_a.py` | `f55168b92b81103a109ed9ae2b423094f0594e5498c93724228d46980005e96d` |
| `tests/test_pb_atlas_numerical_contraction_stage_a_hostile_audit.py` | `c8ccd89bb5ba115d175c7af07d3d89a516b20f87d3a81805cea7bcb0c7135d37` |

## 1. The Schur rule is the exact contraction rule

Let (S=(w_{i_1},\ldots,w_{i_s})) be the original strict classes of curves
already contracted, put (A=S^tGS), and let (c) be another supplied wall.
The pullback of the current image of (c) is its orthogonal projection away
from the exceptional lattice:

\[
       e_c=c-SA^{-1}S^tGc.                         \tag{1}
\]

Consequently

\[
 e_c^2=c^2-(c^tGS)A^{-1}(S^tGc),
 \qquad K_{\rm current}.c_{\rm current}=K.e_c.    \tag{2}
\]

The two executable gates (e_c^2=-1) and (K.e_c=-1) therefore say exactly
that the current irreducible image is a numerical smooth rational
minus-one curve.  The second equality is not redundant on arbitrary arrays:
it excludes a singular or positive-genus numerical image.

Inductively (A) is negative definite with determinant ((-1)^s).  Appending
an eligible curve has Schur complement (-1), so the new Gram remains
negative definite and has determinant ((-1)^{s+1}).  Hence its inverse is
integral.  This proves both termination and the exact integral update used
by the core.  Conversely, every supplied-wall sequential contraction has an
eligible first wall and an eligible next wall after each prefix, so induction
places its mask in the traversal.  Boundary indices are excluded at every
step.

The hostile traversal recomputes each inverse directly, tests equations
(1)--(2), checks every new determinant and all negative-definite leading
minors, and tests **every** remaining allowed wall.  It does not call the
production stabilizer-orbit representative routine.

If the supplied array is the complete negative-curve array of an actual
surface, this also captures every geometric frame: the strict transform on
the original surface of a later minus-one curve has square
(-1-\sum m_j^2<0), so it already occurs in that complete array.  Completeness
and effectivity are hypotheses here; Stage A does not infer them.

## 2. The quotient traversal is complete

The marked permutation group preserves the full wall Gram, all canonical
pairings, and the boundary-weight function.  It therefore preserves both
eligibility equations.  For a reachable representative (M), its stabilizer
partitions the remaining candidates.  Testing one member of each stabilizer
orbit tests every extension of (M).  Replacing an eligible child by its
full-group canonical representative loses no state because the entire child
orbit is reachable.  Induction on depth proves the quotient traversal.

The hostile implementation uses a deliberately different traversal: at each
representative it tests every remaining wall, canonicalizes only the
eligible children, and sums their independently recomputed orbit sizes.  It
reproduces every frozen orbit count and labelled count in all four cases,
including the `H17/48` peak of `311238` labelled masks.

### Repaired supplied-group gate

Before the repair, `validate_automorphisms` accepted, for example, the
identity together with an order-three marked-Gram automorphism while omitting
its square.  Both rows were individually valid, but that two-row set was not
a group; a candidate “orbit” could then escape the alleged stabilizer
complement.  The validator now:

1. rejects duplicate rows and malformed permutations;
2. requires the identity;
3. greedily generates the subgroup of the supplied permutations;
4. rejects any generated product outside the supplied set; and
5. requires final equality between the generated subgroup and supplied set.

This also proves inverse closure, since a finite composition-closed set of
permutations containing the identity is a group.  The algorithm costs
(O(|G|r g)), where (g) is the small number of discovered generators,
rather than constructing the quadratic multiplication table.  Hostile tests
exercise both a missing product and a duplicate identity.  Every frozen
group of orders `24`, `8`, `120`, and `1152` passes unchanged.

## 3. Nonsaturated wall spans are handled correctly

When the wall columns generate the Picard lattice integrally, a marked-Gram
permutation defines an integral lattice map: Gram preservation respects all
relations, the inverse permutation supplies its inverse, and canonical
pairings force (K) to be fixed.  This justifies the saturated shortcut.

For a proper full-rank wall sublattice, the same conclusion is false.  The
core selects a rational wall basis, solves for the unique rational lift, and
retains the permutation exactly when the lift is integral unimodular, maps
every wall correctly, preserves (G), and fixes (K).  Integral lifts are
closed under products and inverses, so the filtered rows remain a group.

The physical `H10/30` eleven-wall subset has full rational span and index
two.  It parses successfully; identity-only and filtered-symmetry traversals
both give the labelled level census

```text
1, 6, 11, 12, 10, 5, 1, 0, 0, 0, 0.
```

The separate synthetic index-two packet has a marked-Gram swap whose unique
Picard lift contains denominator two.  An independent rational solve
recovers that half-integral lift, and the core discards the swap while
retaining the identity.  Nonsaturation therefore reduces compression; it
does not reject the branch or merge inequivalent frames.

## 4. Endpoint reconstruction is integral and complete

For a terminal mask, the contracted Gram is unimodular and the orthogonal
projector

\[
       P_S=1-S(S^tGS)^{-1}S^tG
\]

is integral.  Its image is exactly (S^\perp\cap L), because the exceptional
lattice is a unimodular direct summand.  Column HNF therefore returns the
whole endpoint lattice, not a finite-index sublattice.

In endpoint rank one, nef orientation together with (H^2=1) and
(K.H=-3) gives the numerical plane basis.  In rank two, the two extreme
rays of all projected wall inequalities are recomputed.  If they can be
ordered as (f,h) with

\[
 f^2=0,quad f.h=1,quad K.f=-2,quad
 h^2=e\ge0,quad K.h=-(e+2),
\]

then (s=h-ef) gives the exact Hirzebruch basis

\[
 (s,f)^tG(s,f)=\begin{pmatrix}-e&1\\1&0\end{pmatrix},
 \qquad K=-2s-(e+2)f.
\]

For each sequential order, equation (1) recovers mutually orthogonal total
exceptional transforms (E_j) with (E_j^2=K.E_j=-1).  The emitted basis
((H,E_j)) or ((s,f,E_j)) is checked to be unimodular, and its Gram and
canonical coordinates are respectively

```text
diag(1,-1,...,-1),             (-3,1,...,1),
diag([[-e,1],[1,0]],-1,...),   (-2,-e-2,1,...,1).
```

The hostile verifier independently reconstructs all of this from the frozen
input arrays for all `282` terminal orbit representatives.  It recomputes
every total exceptional vector, basis determinant, Gram, canonical vector,
wall coordinate table, boundary row, endpoint image, and nef ray.  It also
expands every representative under the frozen group and matches all `4865`
labelled terminal masks exactly.

## 5. Four exact calibration censuses

The independent all-candidate traversal gives:

| packet | group | ruled labelled endpoints | plane endpoints |
|---|---:|---|---:|
| `H10/30` | 24 | `F0:837, F1:2064, F2:390` | 573 |
| `H11/26` | 8 | `F0:86, F1:212, F2:54, F3:24, F4:12` | 52 |
| `H14/60` | 120 | `F0:105, F1:310, F2:120, F3:20` | 5 |
| `H17/48` | 1152 | `F0:1` | 0 |

It reproduces the complete labelled depth trails printed in the theorem,
not only the four terminal histograms.  The aggregate remains

```text
F0:1029, F1:2586, F2:564, F3:44, F4:12, P2:630.
```

Unimodular basis mutation and full wall reversal also replay.  The input and
output schemas validate, the fixture extractor reproduces byte-for-byte,
and the regenerated archive pins the repaired core hash.

## 6. Exact scope

This PASS is an all-contraction theorem **relative to the supplied numerical
wall array**.  It neither proves that the walls are effective irreducible
curves nor that the array is the complete Mori/facet array.  It does not
construct a surface, generate Stage-B incidence equations, prove a nonempty
regular open, exclude hidden walls, return the same Coxeter polytope, or give
an ABSW-type dimension theorem.  No Nikulin--Vinberg classification or
Prasad machinery enters the proof or certifier.

## Replay

```bash
/Users/sasha/.venvs/global/bin/python \
  scripts/freeze_pb_atlas_stage_a_fixtures.py --check

/Users/sasha/.venvs/global/bin/python \
  scripts/certify_pb_atlas_numerical_stage_a.py \
  --check-output data/enumerations/pb_atlas_numerical_contraction_stage_a.json

/Users/sasha/.venvs/global/bin/python -m pytest -q \
  tests/test_pb_atlas_numerical_contraction_stage_a.py \
  tests/test_pb_atlas_numerical_contraction_stage_a_hostile_audit.py
```

The combined focused suite has 13 passing tests: eight owner tests and five
independent hostile tests.  The full replay, including all basis/wall
mutations and the separate no-stabilizer traversal, passes.
