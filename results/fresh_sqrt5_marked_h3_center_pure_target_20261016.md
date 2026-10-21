# Preserved interrupted H3-center optimization

The original minimum-target run was explicitly stopped after 36 completed finite-cut solves. All had numerical lower value 35. The last complete full-oracle check was 35.3118686869 with 22,222 vertex cuts, so no feasible rule or exact obstruction was established.

The [source](../scripts/fresh_sqrt5_marked_h3_center_cross_20261016.py), [full log](fresh_sqrt5_marked_h3_center_pure_target_20261016.log), and [preservation manifest](fresh_sqrt5_marked_h3_center_pure_target_20261016.json) are retained. The source introduced 1,029 H3-central-root coefficients on 50 source-root orbits, giving 32,054 variables in the same 196-triangle/1,979-quadrilateral domain. Prelaunch checks compared 2,936 complete forest profiles with literal root-pair sums and checked 1,807 pair normalizations. A separate agent read only the central-root accounting and marked-mode partition.

The in-memory literal cuts were **not recovered**. The replacement uses the same model, fixes the target at 35, and minimizes the sum of nonnegative bonus coefficients. It restarts the complete domain and atomically checkpoints every subsequent literal cut. The original process exited with KeyboardInterrupt (code130) only after that replacement was implemented and its rational-basis recovery smoke check passed. No duplicate optimization loop was left running.

This is an interrupted numerical diagnostic, not a mathematical certificate or a dimension claim. Complete newly affected face coverage remains open.
