# Terminal compact `H9/14` CEGAR audit: Record 0, `a=1`, source 3

Status date: 2026-08-08.

## Exact result

The independently rebuilt terminal CP-SAT model is **infeasible**.  Therefore
there is no compact hyperbolic Coxeter polytope in `H9` with 14 facets that
simultaneously has

- missing-edge graph Record 0, graph6 `M?????????????o??`;
- Dehn--Sommerville schedule `a=1`, meaning one missing order-four support;
- the fixed source-orbit-3 representative `{2,3,4,13}` for that support; and
- the required compact incidence, Lanner, and elliptic vertex-label
  conditions.

This is a terminal theorem for one fixed source, not an exploratory timeout
or a count of failed samples.

## The boundary of the theorem

The quantifier must not be enlarged:

| question | status |
|---|---|
| this fixed Record-0 / `a=1` / source-3 compact branch | **impossible** |
| the other three Record-0 / `a=1` source orbits | not decided here |
| other schedules or missing-edge records | not decided here |
| all compact `H9/14` Coxeter polytopes | not decided here |
| finite-volume noncompact `H9/14` polytopes | **not checked** |
| global Gram realization after a surviving label assignment | unnecessary on this branch because no label/incidence target survives |

In particular, noncompact finite-volume polytopes may have ideal vertices and
affine vertex diagrams.  The present compact elliptic-vertex contradiction
cannot be repurposed as a noncompact volume theorem.

## No-trust reconstruction

The audit script does not import the exploratory CEGAR driver, its dual model
builder, or its core extractor.  It uses only the frozen checkpoint, the
frozen global support schedule, standard Python, and OR-Tools.

Starting with the graph6 string, it independently reconstructs the two
missing edges `{0,13}` and `{1,13}`, the fixed A4, its eight forced A5 shadow
members, the schedule of 32 missing independent five-sets, and 178 scheduled
rank-nine vertices.  It then rebuilds the complete necessary compact dual
model:

| model component | count |
|---|---:|
| candidate facet complements | 824 |
| face-witness constraints | 1,324 |
| all-order Lanner-closure constraints | 7,051 |
| ridge degree-zero-or-two constraints | 1,547 |

These counts agree exactly with the checkpoint.  Even if this incidence
model is viewed merely as a necessary-condition relaxation, terminal
infeasibility is decisive: a genuine compact polytope would give a model of
the relaxation.

## Why each local core is a theorem

All 41 stored cores were replayed independently.  For every core, the audit
does the following.

1. It reconstructs every labelled placement of the classical Lanner
   supports: all 108 placements of the nine order-four types and all 420
   placements of the five order-five types.
2. It merges every mutually compatible placement, including all overlaps
   between supports.
3. It checks every stored face restriction against a separate exact
   classifier for finite Coxeter diagrams of types `A`, `B`, `D`, `E`, `F`,
   `H`, and `I2(3,4,5)`.
4. At every stage it reproduces the checkpoint's compatible-extension,
   restriction, and passing-state counts exactly.
5. The terminal sparse-label state count is zero for all 41 cores.

The enumeration is complete for the target alphabet `{2,3,4,5}`.  An
unassigned meeting pair is initially given label two.  In the Coxeter Gram
matrix this minimizes the corresponding nonnegative cosine adjacency entry;
increasing it cannot decrease the spectral radius and therefore cannot turn
a nonelliptic restriction elliptic.  Thus a failed minimal completion cannot
be rescued by assigning a larger label.

The audit also checks the literal cut semantics.  A face variable is true if
and only if at least one selected rank-nine vertex contains the stored
restriction.  Projecting those cut restrictions to the final support union
gives exactly the independently replayed elliptic restrictions.  Hence the
conjunction forbidden by each clause is genuinely incompatible with compact
Coxeter labels.

## Why the full-source cuts are lossless

For the fixed graph and A4, the complete source stabilizer is

\[
 S_{\{2,3,4\}}\times S_{\{5,\ldots,12\}}\times
 S_{\{0,1\}} \cong S_3\times S_8\times C_2,
\]

of order 483,840.  Adjacent transpositions in the first two blocks and the
swap `0 <-> 1` give ten generators.  The audit verifies that every generator
preserves the missing graph, fixed A4, missing-five ground, candidate facet
ground, and all model constraints.

For every local core, generator-orbit BFS rebuilds its complete orbit.  Every
orbit size and orbit-stabilizer quotient agrees with the checkpoint.  The
orbit-size histogram is

```text
84 x 1, 336 x 2, 1008 x 2, 2016 x 16, 2520 x 2,
4032 x 1, 5040 x 2, 10080 x 5, 20160 x 8, 80640 x 2.
```

The 41 orbits supply 427,140 exact image clauses.  The cut-arity histogram is

```text
5:   5,040
6:  56,784
7: 261,492
8:  55,440
9:  48,384
```

Every image of a valid compact-label contradiction is again a valid
contradiction.  No target is removed on heuristic similarity.

The base model uses one lexicographic inequality for each source generator.
This slice is lossless: the global lexicographic minimum of any finite source
orbit satisfies all generator inequalities.  Therefore infeasibility of the
sliced, fully cut model closes the whole fixed source, not merely its chosen
label order.

## Terminal model and solve

After installing all cores, the independently constructed terminal model is

| quantity | exact value |
|---|---:|
| Boolean/integer variables | 10,602 |
| constraints | 458,793 |
| Boolean-or constraints | 437,124 |
| linear constraints | 20,122 |
| table constraints | 1,547 |
| serialized model size | 36,743,797 bytes |
| binary model SHA-256 | `7485e42496c5f9f8f03cedde468a9f7180b98c9ee92ab6f86a23e68a0b735ec3` |

OR-Tools `9.15.6755` solved this exact model with one worker, random seed zero,
and no time limit:

```text
status                  INFEASIBLE
solver wall time        565.535626 s
total independent replay 640.843332 s
branches                2,783,503
conflicts                 704,576
Boolean propagations  932,148,310
integer propagations   54,841,428
deterministic time        751.355952
```

One exact no-timeout terminal solve was run.  The routine regression suite
does not waste another nine and a half minutes: it independently rebuilds all
41 label proofs, all full source orbits, all 427,140 clauses, and the exact
terminal proto hash, then compares those structural results with the frozen
infeasibility run.

## Pinned inputs and reproducibility

- checkpoint SHA-256:
  `f0b38d563b0b532f39af73203e40050cfe65314f3db50894cac84e4a7f27ef7a`
- support-schedule SHA-256:
  `8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d`
- [independent verifier](../scripts/audit_h9_14_record0_a1_source3_terminal_cegar.py)
- [frozen certificate](../data/enumerations/h9_14_record0_a1_source3_terminal_cegar_audit.json)
- [regression tests](../tests/test_h9_14_record0_a1_source3_terminal_cegar_audit.py)

Full theorem replay:

```bash
PYTHONPATH=src:. python3 -u \
  scripts/audit_h9_14_record0_a1_source3_terminal_cegar.py \
  --check-archive
```

Structural regression without a redundant terminal solve:

```bash
PYTHONPATH=src:. pytest -q \
  tests/test_h9_14_record0_a1_source3_terminal_cegar_audit.py
```

The structural suite passes `4/4` tests in 72.39 seconds.  Frozen artifact
hashes are

```text
independent verifier 6540c3c39afb6de6a710f50348c2ff0bc2d4a42cf2bda72d4ce96c37aea219cf
certificate          afe350d601d76707026475c328479c99c4080809e4665a30b585b696fa61f1c4
regression tests     ad930228c382eadb122dcfdbfc4fdf3e9308b29ec25c53a512bdc776e21f8c0e
```
