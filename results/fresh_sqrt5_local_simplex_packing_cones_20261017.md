# Fourteen exact local simplex-packing optima

The complete local disjoint-block packing cone gives **no additional mass
increase** on the fourteen examined vertex forests, after the frozen subset
and intersecting-orbit clique corrections. This conclusion has matching exact
charge and packing-mixture certificates, rather than only unsuccessful cut
searches.

The repaired partial-model dual therefore remains

\[
 {2790968068837079040\over217745752925063023}
 \simeq12.81755456235>12.
\]

The first ten forests were selected by the largest vertex-mass increases in
the [complete subset scan](fresh_sqrt5_all_subset_port_packing_20261017.md).
The next four were selected because the subsequent
[intersecting-orbit clique scan](fresh_sqrt5_intersecting_block_cliques_20261017.md)
found actual violated inequalities there. They are a justified diagnostic
sample, not a complete examination of the other 726 forests.

## What is optimized locally

Fix one full rank-13 finite forest. Each orbit j consists of literal port
blocks belonging to maximal simplex corners. Let w_j be its total required
corner mass in the saved dual. By the
[disjoint-port theorem](fresh_simplex_disjoint_ports_20261017.md), the blocks
present at an actual vertex are pairwise disjoint.

Let P range over every packing of the selected literal blocks, and let a(P)_j
count its blocks in orbit j. Arbitrary nonnegative orbit charges q_j have
vertex cost

\[
 c(q)=\max_P\sum_j a(P)_j q_j.
\]

The local necessary vertex mass is the exact finite LP

\[
 \max_{q\ge0,\ c(q)\le1}\sum_j w_jq_j
 =\min_{\lambda_P\ge0}
 \left\{\sum_P\lambda_P:
             \sum_P\lambda_Pa(P)_j\ge w_j\text{ for every }j\right\}.
\]

The certificate supplies both sides with equal rational values. A charge is
checked against every packing count vector; its maximum cost is exactly one.
Every positive mixture coefficient comes with explicit disjoint literal
blocks. Their exact weighted orbit counts cover all required marginals and
their total mass equals the charge's exact load.

Orbit counts suffice here. Averaging a literal packing over the finite-forest
automorphism group gives uniform marginals within each complete orbit. The
same averaging can be applied simultaneously to every orbit. Conversely,
averaging any literal charge preserves its uniform-marginal objective and
cannot increase its maximum packing cost. Thus the certificate concerns the
whole nonnegative packing-charge cone on these selected blocks, not only
subset-floor charges or clique charges.

Unknown blocks receive zero charge and can be omitted from a packing. The
relaxation allows combinatorial packings that might fail additional Coxeter
geometry. It does **not** assert that the packing mixtures are actual ambient
vertex stars, or that independently feasible local mixtures fit together
globally.

## Complete finite enumeration and exact results

For each available-port mask, take its least port p. A packing either leaves p
unused, or contains exactly one selected block through p. Removing p or that
block gives a complete recursion on smaller masks. The implementation stores
every achievable orbit-count vector; no heuristic pruning is used.

Across the fourteen forests this produces 10,886 distinct count vectors.
The largest individual list has 3,448 vectors. The exact primal certificates
use 140 positive packing mixtures in total. All numerical proposals were
recovered rationally and checked independently of solver tolerances.

For the first ten forests, their full local optimum equals the preceding
all-subset repaired mass. Their indices are
`536, 289, 503, 2989, 3113, 397, 986, 3317, 910, 3120`.
The saved JSON records their full finite forests, every block orbit, literal
mask and exact optimum.

The four additional cases establish that the newly found clique corrections
already reach the full local packing optimum there:

| Vertex index | Full finite forest | Exact local optimum |
| --- | --- | --- |
| 879 | `3A1+A2+D4+H4` | 2309207 / 5000000000 |
| 904 | `3A1+A2+2H3+I2(5)` | 32925897 / 10000000000 |
| 1451 | `2A1+A2+A3+A4+I2(5)` | 2690016221 / 400000000000 |
| 2948 | `A1+2A3+A4+I2(5)` | 15118781 / 2500000000 |

Their respective increases over the all-subset masses are
`3047/30000000000`, `42409/1000000000`,
`532060553/13200000000000`, and `45697/100000000000`.
There is no further increase over the clique-repaired masses.

## Interface, verification, and scope

- [Exact source](../scripts/fresh_sqrt5_local_simplex_packing_cones_20261017.py)
- [Saved local primal/dual certificates and charge columns](fresh_sqrt5_local_simplex_packing_cones_20261017.json)
- [Pointed-base redundancy certificate](fresh_sqrt5_pointed_simplex_base_capacity_20261017.md)

`local_cones` contains the fourteen complete local certificates.
`exact_charge_columns` gives a cost-one column per forest, with rational
coefficients on the frozen positive simplex rows. Quad coefficients are zero.
The full saved alpha is unchanged; the final mu is exactly the preceding
clique-repaired vector. Thus all old positive-cost constraints, signed
transfers, actual provenance assignments, and pointed-base constraints remain
valid. The unchanged exact bound above is a partial-model obstruction only.

Generation uses fourteen small local solves, each with a twenty-second cap.
There is no global optimization. Run
`python3 scripts/fresh_sqrt5_local_simplex_packing_cones_20261017.py --verify`
to reconstruct the complete packing lists and verify both exact certificates
without invoking an optimizer. The `--refresh` option likewise reuses the
saved exact certificates and performs no solve.

**FINAL/FROZEN.** Generation, a source-metadata refresh from the exact saved
certificates, and a separate fresh saved replay passed. Direct source hashes,
Python syntax, and local report links passed. These are internal checks; no
complete geometric face coverage, global realizability, external refereeing,
or new dimension bound is claimed.
