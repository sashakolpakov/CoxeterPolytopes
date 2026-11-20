# Linear programming for dimension bounds

The completed LP arguments recover **57** for finite volume, **17** for
compact polytopes and **15** for compact arithmetic polytopes. They combine
proved local inequalities on exhaustive configuration classes. Deriving those
local theorems from elementary Gram conditions, or improving the bounds,
remains separate work.

Read these in order:

1. [Short reading guide](LP_READING_GUIDE.md): the variables, primal and dual
   problems, diagram marks, and dimension monotonicity.
2. [The uniform LP formulation](UNIFORM_LP_PARADIGM.md): the common
   certificate argument and the different inputs in each of the three cases.
3. [Completed recoveries and execution evidence](existing_proof_lp_recovery_results_20261118.md):
   the exact results and their limitations.
4. [The 36-page manuscript](../../arxiv/linear_programming_dimension_bounds/main.pdf)
   and [submission guide](../../arxiv/linear_programming_dimension_bounds/README.md):
   the detailed proofs, source files and typesetting records.

The formulation notes give the individual models:

- [Common energy framework](existing_proof_energy_lp_master_20261118.md).
- [Compact bound 17](compact17_complete_energy_lp_recovery_20261118.md).
- [Compact arithmetic bound 15](compact_arithmetic15_global_energy_lp_recovery_20261118.md).
- [Execution contract](existing_proof_energy_lp_execution_contract_20261118.md).
- [Manuscript language revision and primary-source reading notes](manuscript_language_revision_20261119.md).

The compact and arithmetic model exports and exact certificates accompany
these notes. The finite-volume record also depends on previously retained
finite data; its complete older matrix and the repository's unuploaded bulk
collection are not included in this publication. A build receipt establishes
typesetting, not a new mathematical replay.

Earlier exploratory LP/SDP work and its unsuccessful runs are retained
locally. They do not replace the completed existing-proof recoveries or
establish stronger bounds. The full working index has been preserved
byte for byte in the local session history; this public index selects the
material needed to read the current manuscript.
