# General endpoint scouts in the mixed H9 chamber

The numerical search produced two exact compatible endpoint walls, certified separately in the [fourteen-wall construction](fresh_h9_mixed_d4_endpoint_pair_20261014.md). That chamber still has infinite volume. The numerical banks described here are proposals, not exact Coxeter certificates or complete nonexistence results.

The [first scout](../scripts/fresh_h9_mixed_d4_endpoint_caps_20261014_scout.py) uses the [complete local profile census](fresh_h9_mixed_d4_endpoint_profiles_20261014.md) on all 52 open directions of the twelve-wall chamber. It tests 1,825 finite profiles and nine labels in each of 72 isolated dihedral tails, totaling 2,473 sampled profiles. Each continuous all-dotted interval contributes one sample. Finite-angle anchor equations use labels 2,3,4,5,6,7,8,9,10,11,12,16,20,24,32,64; all remaining wall angles are screened numerically against labels 2 through 256 or a strict dotted inequality.

The [saved first bank](fresh_h9_mixed_d4_endpoint_caps_20261014_scout.json) has 313,664 real anchor branches and 35,486 distinct accepted numerical proposals. Of these, 4,279 cut the original strict escape. Proposals 98 and 324 give the two independently verified exact walls. Their exact validity does not depend on this numerical search.

The [next scout](../scripts/fresh_h9_mixed_d4_endpoint_pair_20261014_next_caps.py) searches the 21 open directions starting at the three newly created finite vertices, using the [new complete local profile census](fresh_h9_mixed_d4_endpoint_pair_20261014_next_profiles.md). These directions have 459 finite profiles and 21 unrestricted isolated dihedral tails. Sampling the same nine tail labels gives 648 profiles, 123,744 real anchor branches and 8,928 distinct accepted numerical proposals in the [second bank](fresh_h9_mixed_d4_endpoint_pair_20261014_next_caps.json).

None of the second bank cuts the fixed rational point whose original-coordinate decimals are recorded literally in its source. The point was chosen to survive the first bank as well. This is a numerical search result only: arbitrary tail labels, unsampled anchor labels, the full continuous interval geometry, affine endpoints, and walls not ending an existing finite-corner edge are not excluded by either scout. Increasing wall counts without removing the common escape would not establish finite volume.

Both scouts pass deterministic saved-result replay. Exact geometry and independent arithmetic validation for the selected pair are recorded in its separate producer and review reports.
