# Complete local simplex packing at the new enlarged-model masses

Fresh exact packing covers for every positive-S vertex forest lower the enlarged-model dual from12.52181135130475… to

\[
\frac{8460681930118372500000000000000}{680097081479055986799187995041}
=12.440403231430313\ldots>12.
\]

This repaired witness satisfies the entire nonnegative local simplex-block charge family at the saved new face masses, as well as every exported column of the [enlarged finite model](fresh_actual_corner_enlarged_model_20261017.md). It remains a partial-model obstruction, not a dimension theorem. Joint simplex/prism geometry and unexpanded Q provenances are separate obligations.

The [source](../scripts/fresh_actual_corner_enlarged_model_20261017_local_cones.py) reconstructs all2,638 positive S rows at the new alpha. It does not reuse preceding fixed-alpha covers. Their dimensions are1:1,173;2:530;3:480;4:455, giving8,131 vertex-corner occurrences. The complete literal block orbits give5,011 marginal coordinates across815 rank13 finite forests.

For each forest, all disjoint literal block packings are enumerated by the inherited exact thirteen-port dynamic program. One bounded local packing LP is used for rational recovery; no second global LP is invoked. The saved [certificate](fresh_actual_corner_enlarged_model_20261017_local_cones.json) contains672 exact local optima and143 exact cover brackets, totaling4,467 positive literal packing-mixture terms. No resource guard or unresolved local solver occurred. The enumeration limits remain explicit: one million cached count-vector entries,30 seconds enumeration and20 seconds local solve per forest; a resource stop would retain a valid singleton cover and mark optimality open rather than invent infeasibility.

The new local separation finds115 exact violated charge columns. Exact covers require215 vertex-allowance increases. Every face mass, including all2,294 Q masses, stays fixed. All430,666 exported model columns are directly checked after repair, including96,176 zero-budget signed columns and52 signed vertex-cost incidence columns.

## Why finite mixtures cover the entire local charge cone

At a finite vertex, the port blocks of distinct maximal simplex faces are disjoint by the frozen [disjoint-port theorem](fresh_simplex_disjoint_ports_20261017.md). Each saved mixture term is an explicitly checked packing of disjoint thirteen-bit masks. Its rational weight is positive, and the weighted orbit counts cover every required new-alpha marginal.

Averaging each **whole packing simultaneously** over the finite forest automorphism group gives uniform marginals within every complete literal block orbit and preserves disjointness. Thus, for any nonnegative literal block charge, its total marginal load is bounded by the mixture mass times its maximum charge on a disjoint packing. The repaired vertex allowance dominates this mixture mass. This proves all such local charge inequalities at once, not merely the115 exported separating charges.

The argument uses upper covers, so an exact local optimum is unnecessary in the143 bracket cases. Block orbits absent from the positive face support have zero required marginal; allowing charges on them cannot invalidate the cover argument. This does not assert that arbitrary abstract disjoint block packings extend to actual ambient geometric configurations.

## Verification and scope

Generation8108 and fresh saved replay81478 completed PASS. The replay reconstructs every literal orbit and packing domain, verifies each saved rational charge and mixture without optimizing, compares the entire saved result, and rechecks the full finite matrix. Source SHA256 is `915aedcdb06e9c9d4540bcfb22790b2af8e6ce7d929dcdc3e84b5ca562ddcd04`; JSON SHA256 is `07ba3dbf957fde2343f15feefcfb62d8aaab5b5e2ffe89f9be8c80f9b09be306`.

Construction-audit's narrow source/accounting review passes the new-alpha feature reconstruction, progress-cache binding, literal mixture checks and complete430,666-column recount, including the52 signed vertex-cost columns. It does not duplicate the arithmetic replay or re-audit the whole inherited solver engine.

This is complete local simplex-packing coverage for the new positive S support. It does not supply geometric representatives for the1,897 unselected merged Q rows identified in the [support report](fresh_actual_corner_enlarged_model_20261017_support.md), and it does not transfer old joint-port scans to changed masses. A new joint scan treats the397 actual Q rows explicitly and uses zero new coefficients on unexpanded classes.
