# Seven selected H9/14 nerves are not compact Coxeter nerves

## Re-audited theorem

The seven explicitly displayed ordinary-polytopal nerves numbered

```text
30, 32, 35, 42, 43, 44, 47
```

are not nerves of compact hyperbolic Coxeter 9-polytopes.

This statement concerns exactly those seven selected nerve completions. A
record number denotes a missing-edge graph shadow, so the theorem does not by
itself exclude another completion of that shadow. It also does not exclude a
noncompact finite-volume chamber with different ideal-face data.

## Scope re-audit

On 2026-08-06 a hostile review provisionally demoted the pairwise
disjoint-Lanner arguments for records `30,42,43,44,47`, on the theory that
outside dotted or Lanner components might mediate between the selected
supports. A second hostile audit rejected that concern.

The decisive matrix is a **principal Gram block** on the union of two
disjoint Lanner supports. Entries incident to nodes outside that union cannot
change the block or its inertia. Therefore outside components cannot repair
two negative directions already present in that block. Likewise, a failed
actual-vertex restriction is a principal subdiagram of the required elliptic
vertex diagram; adding or strengthening non-right entries cannot make it
elliptic.

The new certificate does not merely restore the former hash. It explicitly
freezes the support-symmetry orders, the complete cross-pair case split, and
the outside-entry principal-block lemma. The separate audit note is
[`h9_14_selected_support_scope_reaudit.md`](h9_14_selected_support_scope_reaudit.md).

## General disjoint-Lanner gate

Let `L` and `R` be disjoint minimal nonfaces. In a compact Coxeter chamber,
their induced diagrams are Lanner. Every pair in `L x R` below is a meeting
pair occurring in an actual rank-nine vertex, so it carries a finite label.
All labelled Lanner placements are covered either literally or modulo a full
support symmetric group verified on the complete maximal-face family.

For a fixed placement pair:

1. If there is no non-right `L x R` edge, then
   `G[L union R] = G[L] direct-sum G[R]` has negative index two.
2. A single edge is first tested at the weakest label `m=3`. If an actual
   vertex is already nonelliptic, Perron--Frobenius monotonicity rules out a
   stronger label and every super-labeling.
3. Every individually admissible edge is tested at `m=4,5`, and every pair
   of individually admissible `m=3` edges is tested. Any larger edge set
   contains one of these monotonicity-minimal bad subsets.
4. Every remaining one-edge principal block has exact negative index two.

For records `43,44,47`, fix any one of the three audited residual A5
supports. The zero-, one-, and at-least-two-edge cases exhaust all entries
between it and the normalized A5. Edges with an endpoint outside that fixed
ten-node block cannot change it. If another residual support overlaps the
fixed one, any edge whose two endpoints lie in the fixed block is already one
of its audited entries. Thus even a fixed support gate is decisive; all three
are retained as redundant checks.

## Exact census

| selected nerve | rank-nine vertices | verified normalization | decisive census | outcome |
|---:|---:|---|---|---|
| 30 | 128 | `S4 x S4`, order 576 | `9*9*16 = 1,296` single-`m=3` trials | zero survivors; orthogonal A4+A4 otherwise |
| 32 | 110 | complete overlapping-A4 join | `48*420*20 = 403,200` trials | zero survivors |
| 35 | 102 | independent common-four-core certificate | `420 -> 1,332 -> 4,620`; 96 rank-five survivors | all killed by affine actual-vertex witnesses |
| 42 | 100 | `S4 x S5`, order 2,880 | `9*5*20 = 900` trials | two isomorphic one-edge exceptions |
| 43 | 90 | full left `S5`, order 120; right placements literal | three gates of `5*420*25 = 52,500` | 300 one-edge exceptions per gate |
| 44 | 90 | same | same | same |
| 47 | 90 | same | same | same |

The 108 labelled order-four Lanner placements split into nine `S4` orbits
of sizes

```text
6, 6, 12, 12, 12, 12, 12, 12, 24.
```

The 420 order-five placements split into five `S5` orbits of sizes

```text
60, 60, 60, 120, 120.
```

These are precisely the canonical placement lists used when a support is
normalized. Other supports are retained literally whenever their full
symmetric group is not a nerve symmetry.

For record 42, the two exceptional `m=3` trials occur in one weighted-block
class. Labels `m=4,5` leave zero survivors. Of 8,550 raw two-edge trials,
monotonicity leaves only the pair of the two exceptional edges, and it fails
an actual vertex. The remaining nine-node block has exact leading-minor signs

```text
++++++-+
```

and inertia `(7,2,0)`.

For each support gate of records `43,44,47`, 300 `m=3` trials survive, while
`m=4,5` and every residual double-edge trial have zero survivors. Up to
weighted-block isomorphism the exceptions form two classes. Both have exact
leading-minor signs

```text
++++-----+
```

and inertia `(8,2,0)`.

All radical signs are proved from rational enclosures for `sqrt(2)`,
`sqrt(5)`, and `sqrt(10)`; no numerical eigenvalue decision is used.

## Records 32 and 35

Record 32 uses a different complete route. The two overlapping A4 supports
have 546 compatible joins and 48 states after all actual-vertex restrictions.
All 420 placements on the disjoint A5 are then retained, and all twenty cross
pairs are tested. The 403,200 minimal trials have zero survivors; the
orthogonal alternative has negative index two.

Record 35 is chained to its independent common-four-core certificate. Its 96
intermediate states form two symmetry orbits, with affine `D4+A1` and affine
`D5` witnesses inside actual vertices. Entrywise monotonicity makes that an
unrestricted compact-label obstruction.

## Conclusion and finite-volume warning

Every labeling of each selected nerve contains either a monotonicity-stable
nonelliptic actual-vertex restriction or a principal Gram block with at least
two negative directions. Cauchy interlacing contradicts the single negative
direction of a hyperbolic Gram matrix.

This is not a finite-volume theorem. A noncompact candidate still requires
the literal Vinberg/Guglielmetti all-edge endpoint test and the complete
connected-parabolic zero-neighbourhood audit. Cusp counts alone do not
suffice.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_h9_14_selected_nerve_coxeter_no_go.py --verify
PYTHONPATH=src:. python3 scripts/independent_h9_14_selected_nerve_coxeter_no_go_audit.py
PYTHONPATH=src:. python3 -m pytest -q tests/test_h9_14_selected_nerve_coxeter_no_go.py
```

The machine-readable certificate is
`data/enumerations/h9_14_selected_nerve_coxeter_no_go.json`.
