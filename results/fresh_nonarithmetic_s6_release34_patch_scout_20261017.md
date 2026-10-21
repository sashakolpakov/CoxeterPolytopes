# Releasing one old contact closes the numerical three-cap equations

The [scout](../scripts/fresh_nonarithmetic_s6_release34_patch_scout_20261017.py)
deletes old walls 34 and 52 from the 63-wall chamber, retains the other 61,
and replaces 52 by a continuously movable normal. Its contacts with
`S6=[0,1,2,3,5,6]` and old 8 stay right, and its contact with old 30 stays
pi/4. The former right contact with deleted 34 is released.

The norm-two completions form a circle in the positive two-dimensional
orthogonal complement of `S6+[8,30]`. Writing its particular vector as p
and equal-length orthogonal residual vectors as u and v gives

`r52(t) = p + ((1-t^2)u + 2tv)/(1+t^2)`.

For each t, C has right contacts with S6 and 8, pi/15 with 44, and pi/4
with replacement 52. A and B are the negative normalized polars of 52 and
44 in `[30,44,replacement52,C]`. Thus A is right to S6,30,44,C and B is
right to S6,30,replacement52,C. Their remaining product is tuned to pi/27.

The [numerical result](fresh_nonarithmetic_s6_release34_patch_scout_20261017.json)
finds one root near `t=-5.1241142736e-8`. The A--B product residual is about
`1.7e-12`; all 250 new contacts pass the numerical gates, comprising 36
prescribed finite contacts and 214 strictly dotted contacts. Eight proposed
finite corners are positive definite and feasible. The old center remains
strict. Replacement 52 has product approximately `-1.36e-6` with deleted
34, so this calculation does not retain that wall.

Generation and a fresh saved replay passed. These are floating-point discovery
results. Exact root isolation, contacts, actual facets and complete local
polygons are established separately by the
[exact construction](fresh_nonarithmetic_s6_release34_exact_20261017.md).
The [boundary comparison](fresh_nonarithmetic_s6_release34_boundary_20261017.md)
also accounts for the boundary lost when old 34 is deleted. This scout does
not establish finite volume or an all-parameter exclusion.
