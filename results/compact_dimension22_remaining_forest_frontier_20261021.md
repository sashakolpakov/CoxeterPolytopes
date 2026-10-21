# Dimension22 after the uniform pressure refinements

21 October 2026. Complete written root composition for separate review.
The new component proofs have full source-bound mathematical reviews.
No scientific, numerical, Lean or TeX job ran. The global compact
ceiling remains22; this note identifies what is still required to
lower it. It is supporting work for the same compact paper.

**Proposition.** If a compact hyperbolic Coxeter22-polytope exists,
then at any vertex farthest from any strict interior point its finite
diagram must be one of

    B22, D22, B21+A1, D21+A1.

In the last two cases the A1 has an isolated short opposite of type I.
The statement is about farthest vertices; it does not assert that every
vertex can be made farthest from an interior point. No arithmeticity
is assumed. These four alternatives have not been excluded here.

## 1. Accepted inputs and decisions in this batch

All paths below are under results/. Source-bound mathematical proof,
byte integrity and machine execution are different checks.

| Input or new theorem | Source and SHA256 | Written review |
| --- | --- | --- |
| Old complete compact22 ceiling and baseline component accounting | `compact_dimension23_exclusion_composition_candidate_20261020.md`, `8a691e4fb5e51f3256a02c68abd70a6bafad3e9dea99ddc45d345585b4d86edd` | Accepted root/global review b314f6a6 and separate review8ccd1c7d; numerical ceiling22 is unchanged. |
| Every A_r,r>=16 has actual energy>9/14 | `compact_large_an_actual_nine_fourteenths_uniform_transfer_candidate_20261021.md`, `cf14513f43a44fa139e6ef0e35ca9ea6f20f0183f397375015ad0ecf5e2d6ab0` | Root author accepts after complete rereading and outside-author review a5aa8022; no source correction. |
| Every A_r,r>=17 has actual energy>31/40 | `compact_large_an_actual_thirty_one_fortieths_uniform_transfer_candidate_20261021.md`, `194caa4fd0d5661fb3b7f579a291cf48da8ea2354391210a9dd1c8210ceaf8b0` | Full root review and outside-author review ac730640; no source correction. |
| Uniform B/D prefix gain G=3969/800000 for M6..16 | `compact_classical_generic_prefix_terminal_pressure_gain_20261021.md`, `194c2b8b4327538d40fc079e137824266274897d773af3b1a598ccde2d27720b` | Full root review099c7aaa and outside-author review; M6 hole wording clarified before acceptance. |
| No A20 with rank-two spectator remainder in dimension22 | `compact_a20_small_debit_coupling_exclusion_20261021.md`, `6ed1ed3c81f5126c6ebf66c07cf8fc3338db3b64352d8985af09db74c6ee9721` | Root author accepts after full rereading and outside-author review942edca2; both I intervals and both7 roles are covered. |

The exact component caps, ordinary rates and defect isolation are the
same accepted inputs as in the baseline. In particular A<=20,B<=22,
D<=23; ordinary components pay more than rank/14; a nonshort A1 pays
more than1/10; at most one short A1 or instead one of the two isolated
rank-two7 defects can occur. The latter has cost>.09. J can be handled
by a one-short-rank discard, or by its stronger combined A1-pair cost>.34.

The generic prefix source is f637d729 and has rank hypotheses B_r>=M+5,
D_r>=M+6. The new gain keeps exactly these hypotheses and has arbitrary
finite-spectator scope. It must not be added independently to another
gain that already fixes the same first-coordinate block.

## 2. A dominant component exists and is unique

At the farthest vertex use the common vector x with ||x||^2<1 and
orthogonal additive component energies. Dominant means A_r,r>=7 or
B_r/D_r,r>=9. With no dominant, finite classification and the ordinary
rates already give a contradiction after either permitted small debit:
21/14>1, or9/100+20/14>1. Every dominant has energy>1/2, so two are
impossible. Only after this uniqueness is proved are all spectators
declared ordinary apart from the named exceptional possibility.

For spectator rank S>=2 the budgets are S/14 with no exception,
(S-1)/14 after a short-rank discard, or9/100+(S-2)/14 for the7 defect.
The last exceeds the short-rank budget by13/700.

## 3. Every A alternative is now impossible

The half bound first gives S<=7, hence r>=15. If r15, the already
accepted A>=8 four-sevenths theorem and S7 give total energy>1 after
either debit. At r16 the new9/14 theorem and S6 similarly give
9/14+5/14=1, strictly exceeded.

At r17 the stronger31/40 theorem and S5 give more than
31/40+4/14>1 after the worst debit. At r18, S4, fully ordinary
spectators are already excluded. After a possible short-A1 discard,
the remaining ordinary rank3 costs more than17/70: an irreducible
rank3 costs>1/4, rank2 plus nonshort A1 costs>1/7+1/10, and three
nonshort A1 cost>3/10. Hence total energy>57/56. With a7 defect,
the remaining ordinary rank2 costs>1/7 and total energy>1411/1400.
This handles J conservatively too, without a false stronger dihedral
floor or overlapping coupling gain.

At r19, S3, the existing6/7 theorem plus the weakest short-discard
budget2/14 already exceeds1 strictly; the7 and ordinary budgets are
stronger. At r20, S2, the new complete coupling source handles every
possibility. Its ordinary and J cases are immediate. Isolated I splits
at height1/5 and uses genuine propagation or its own positive energy.
Either7 defect forces e_A20<.91 while its actual direction floor and
two separated strong tails force e_A20>.9103689. This last step uses
the strengthened near4 estimate, not a conclusion valid only under
the old6/7 contrary target. Higher A ranks are capped. Thus no A
dominant survives in dimension22.

## 4. The exact remaining B/D ranks

Their9/14 rate gives S<=5 and therefore r>=17. For r17,18,19 the
common D-admissible parameters M=r-6 are11,12,13. The corresponding
old component floors .8472225,.8689525,.88864 plus the weakest
spectator budgets4/14,3/14,2/14 already exceed1. The B requirement
is weaker; the new positive gain is not even needed in these rows.

At r20, M14 gives the improved floor

    energy>.90715+G=.91211125.

An ordinary rank2 remainder pays>1/7; a short A1 plus its other
nonshort A1 pays>1/10 even after discarding the short rank; and J
pays>.34. Each exceeds the remaining budget. The only former coarse
exception was the isolated rank-two7 block. Its>.09 now gives total
>1.00211125. Thus r20 is excluded for both B and D.

At r21, parameter M15 applies to both types and gives energy
>.9240775+G=.92903875. The single spectator cannot be J or a
rank-two defect. A nonshort A1's>1/10 already contradicts the norm.
Only isolated short I remains, giving B21+A1 or D21+A1 as stated.

At r22 there is no spectator: whole B22 and whole D22 remain. A
higher rank cannot fit dimension22. This proves the proposition.
It does not assume any arbitrary metric face is Coxeter, transfer a
farthest estimate to a new nonfarthest vertex, or borrow energy twice.

## 5. Construction progress and verification boundary

The separate literal host L=<-phi>+H3+A6 is now excluded at the full
internally reviewed written-proof boundary, including its declared
inherited exact first-corner census. The accepted period source is
`compact_h9_h3_a6_end_edge_period_obstruction_candidate_20261021.md`,
SHA256 `b8985c4b7445f8e21a27775795491a740ba126850f7c7675233e938e104621d3`,
with full root review166350c4 and outside-author reviewdd7669db.
Assumed compactness supplies an actual first end-edge facet; every
possible facet profile either is excluded or forces an integral
full-chamber stabilizer with trace>10. This is a period obstruction,
not a no-finite-endpoint assertion. The explicit standard anisotropic
arithmetic-lattice input then gives noncofiniteness. Other lattices,
the time coefficient2phi host and commensurability classes are not
excluded. No dimension9 example has been obtained.

The H30 positive-height spatial census is source-reviewed and inherited,
not independently rerun in full in this batch. All new statements above
have written exact identities and inequalities; there is no new finite
enumeration, script output, numerical optimization or Lean execution.
The six existing source/PDF pairs and their selected formalization
snapshots remain unchanged. In particular the compact manuscript and
compiled accounting kernel still state23; written22 and this necessary
dimension22 frontier are subsequent supporting research for that one
paper. Lean remains the final step after the new geometry is settled.

The next mathematical task is precisely the four B/D alternatives above,
or stronger arithmetic first-type obstructions, or a genuinely new compact
dimension9 construction. Another global ceiling is not established here.
