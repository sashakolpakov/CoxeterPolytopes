# No current polar wall cuts the48-wall escape

The complete48-wall cone has2,109 hyperideal rays. Their proposed polar hyperplanes split exactly into2,046 whose normals pair strictly negatively with the current escape and63 whose normals pair positively but fail a Coxeter-angle condition. None remains unresolved.

For each of the63 positive directions the certificate exhibits one old wall whose squared normalized product lies strictly between consecutive values `4 cos²(pi/m)`. The gaps are61 between labels2/3, one between3/4, and one between4/5. These are all-label exclusions: the thresholds are strictly increasing for integer `m>=2`; no angle alphabet or cutoff is assumed. Machin's formula with alternating-series bounds encloses pi, and cosine Taylor remainder bounds give outward rational threshold intervals.

Every cone ray is reconstructed through its exact positive clipping expressions. Outward256-bit rational intervals certify all signs and gap inequalities; no exact-field fallback was required in this run. The input is the [complete48-wall cone](fresh_nonarithmetic_joint_regions_20261016.md), and the escape is the exact common witness from its batch certificate. The separate46-proposal numerical screen is not a premise.

This rules out escape-cutting polars of this fixed cone only. Nonpolar walls, changed retained walls, and polars of later chambers remain open. In fact the separate [nonpolar construction](fresh_nonarithmetic_common_escape_20261016.md) finds a compatible cutter.

Generation and a fresh saved replay both passed. [Exact data](fresh_nonarithmetic_common_escape_20261016_polar_gate.json) and [replay source](../scripts/fresh_nonarithmetic_common_escape_20261016_polar_gate.py) bind their inputs by SHA256. This is a primary certificate; no separate independent replay of this gate is claimed here.
