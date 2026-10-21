# Conditioning the three rescalings on the obstructing context still fails

This diagnostic varies the same three old bonus families only when the **actual** hub counts are (A1,I2(5),H3)=(7,2,0). All other bonuses remain fixed. The marked forest 7A1+2I2(5)+2A5 still has exact minimum

54699137087747/1275785690000 ≈ **42.87486332**, above the target 35.

All 1,979 inherited baseline-zero quadrilateral rows are imposed in all 71,036 common spectator contexts. The exact multipliers are 18936/18937 for A5 I2(5)-cross bonuses, 9401/10749 for I2(5) A1-cross bonuses, and 13473/13474 for A5 connected bonuses.

Three literal prior-catalogue quadrilateral rows give matching independent lower bounds: row10 with external counts (5,1,0), row80 with (6,1,0), and row94 with (6,1,0). Their nonnegative rational dual multipliers are recorded in the certificate. The primal satisfies all 4,256 distinct affine profiles, and the dual value equals its objective exactly.

This closes the simple response of changing only the offending context. A useful next proposal must reallocate payments to other corner orbits of the constraining quadrilaterals, change additional bonuses, or change the discount domain. It does not exclude a general compensated rule.

The literal alternative corners make that next step concrete:

- Row10 couples A1×I2(5) in context (7,2,0) to a connected B6 corner in (5,1,0). Its other corners, A1×A4 in (7,1,0) and I2(5)×B4 in (5,2,0), currently receive zero bonus.
- Row80 couples I2(5)×A5 in (7,2,0) to a connected E7 corner in (7,1,0). Its I2(5)×E6 corner in (6,2,0) receives zero. Its fourth corner lies between two separate I2(5) components in (6,3,0), a pair family absent from the inherited hierarchical hub rule.
- Row94 couples connected A5 in (7,2,0) to I2(5)×H3 in (7,2,1). Its H3×A3 corner in (8,1,1) and I2(5)×A4 corner in (6,3,0) currently receive zero.

These are possible reallocation directions, not feasible-rule claims. In particular, introducing a payment between two equal hub components requires its own multiplicity accounting and complete face checks.

The [source](../scripts/fresh_sqrt5_marked_old_bonus_context_scaling_20261016.py) generates and freshly replays the [exact certificate](fresh_sqrt5_marked_old_bonus_context_scaling_20261016.json). See the [global three-family diagnostic](fresh_sqrt5_marked_old_bonus_scaling_20261016.md) and the [complete marked-vertex retained-rule model](fresh_sqrt5_marked_compensation_proposal_20261016.md) for the surrounding scope.

**Scope:** this is an exact necessary-obligation LP for one forest and selected rescaling directions. It has no complete vertex-budget or new affected-face coverage claim, and it changes no dimension bound.
