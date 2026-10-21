# Complete local simplex covers at the new mixed-raw39 dual

The newly optimized finite-model witness remains above the target12 after
checking the complete represented local simplex-packing family. Its exact
repaired ratio is

    4499148564650545000000000000000 / 361055510519686936816210386119
      = 12.46109928685114...

The input [finite model](fresh_mixed_raw39_full_model_20261017.md) had exact
dual ratio12.48512744950078. This is a fixed-witness repair and method
obstruction within the stated domain, not a new global optimization or
dimension theorem. The [source](../scripts/fresh_mixed_raw39_full_model_20261017_local_cones.py)
and [data](fresh_mixed_raw39_full_model_20261017_local_cones.json) are fresh
successors; no frozen predecessor was changed.

## New masses and complete represented features

The new witness has5059 positive face masses:2886 maximal-simplex or edge
rows and2173 quadrilateral rows. Compared with the old4661-face witness,
1341 faces are newly positive,943 disappear and all3718 common positive
faces have changed mass. There are640 newly represented positive simplex
rows. No old fixed-mass packing certificate is assumed valid.

The complete model is rebuilt and its prepared hash checked against the
frozen input. Every positive simplex row is reconstructed, including
isolated maximal edges. Its actual vertex port blocks are expanded through
every weighted-forest automorphism orbit. This gives804 forests,5373 block
orbits and8863 actual corner occurrences. The simplex dimensions1,2,3,4
have respectively1243,662,514,467 positive rows.

For each forest the complete disjoint-block packing domain is enumerated
with the inherited exact dynamic program. Only bounded per-forest linear
programs are used; no global optimizer is invoked. The guards are one
million count-vector entries,30 seconds of enumeration and20 seconds for
each local solver. No guard remains open in this result.

## Certificates and the full local charge family

All804 forests have exact rational covering mixtures. There are700 exact
local optima and104 certified cover brackets, with4658 positive packing
terms containing20074 literal blocks. Each saved term gives explicit
disjoint13-port bit masks.
The verifier checks disjointness, membership in the complete literal orbit,
and coverage of every newly reconstructed orbit marginal.

The disjoint maximal-simplex port theorem makes every actual vertex
configuration a permitted packing. Averaging an entire saved packing under
one common forest automorphism gives uniform literal marginals without
breaking disjointness. Therefore each exact covering mixture bounds the
load of every nonnegative orbit-local simplex charge by its mixture mass
times its maximum packing cost. This is coverage of the whole such local
charge cone on the represented marginals, not only the finite optimized
list of charges. Matching exact charges establish the700 local optima;
the104 brackets still provide valid universal upper covers.

The checker exports78 exact charge columns whose loads exceed the original
vertex allowances. The resulting covers raise164 allowances. The largest
increase is at forest `2A1+2A3+H3+I2(5)` (vertex1665), by
`20486231/80000000000`. All face masses remain fixed.

Every one of the433753 inherited model columns is directly checked after
repair. This includes96176 zero-budget signed transfers and52 signed
incidence columns with a positive vertex cost. Their negative simplex
coefficients are retained. Increasing only vertex allowances preserves
the earlier exported constraints and all unchanged face obligations.

## Verification and remaining coverage

Primary generation (session13551) and a fresh saved replay (session81794)
both exit zero with PASS. The replay reconstructs the domain and verifies
the rational certificates without calling any optimizer. Progress reuse in generation
is keyed to both the exact model hash and the new face-mass vector, so
certificates from different masses cannot be substituted.

This result settles the represented local simplex marginals at this
witness. It does not reconstruct every merged quadrilateral provenance or
all implicit joint simplex/prism constraints at the new masses. The
selected3/4/5-port joint families require their own successor scan.
Broader face, forest and ambient arithmetic compatibility remain separate
obligations; compact-arithmetic dimension21 is unchanged.
