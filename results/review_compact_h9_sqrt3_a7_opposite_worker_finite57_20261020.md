# Independent source review: opposite A7 interval

20 October 2026. **PASS before execution**, for the stated first-edge
computation and its written lattice/chamber inputs. No worker import,
enumeration, numerical job, compilation or typesetting was run here.

Read in full:

* `results/compact_h9_sqrt3_scalar_facet_next_edge_20261020.md`, SHA256
  `f6be0faf51bdc8188a3f8363fe2c5e3ee35e16c356ba49dd1546574b5195567d`;
* `results/compact_h9_sqrt3_a7_opposite_interval_20261020.md`, SHA256
  `6e287ec6d16f9c1700ac717ce605a3ea956afefce1b0768182d711b2d87dcd58`;
* `scripts/check_compact_h9_sqrt3_a7_opposite_edge_20261020.py`, SHA256
  `4ac5bbe521fe0ec6b46eb0c647e2b204c40a78b75e1b4c271029503ef584534a`.

## Lattice and initial corner

The scalar-facet argument explicitly excludes the other norm6 class and all
oblique intersections with the selected norm6 wall. The orthogonal lattice
is therefore the full intrinsic reflection lattice M, rather than a partial
embedded root subsystem. Its integral change of variables and determinant
1-s check. The mod2 argument excludes all unit norms, leaving2 and2epsilon
as the necessary primitive-reflector representatives.

The worker's ambient map has time multiplier423+244s and spatial shift
185+107s. I recomputed

    2(185+107s)(4+s)+9(185+107s)^2-s(423+244s)^2=-2-s,

so it preserves the asserted form on M. The literal orthogonality check is
against rho=(1669+964s;(732+423s)^9), with the correct signs.

Independently, the zero-height chamber of M is exactly A8. If
z_i=A_i+B_i*s and t=0, norm2 gives
`sum A_i^2+3 sum B_i^2=2`, forcing all B_i=0 and roots e_i-e_j. Norm2epsilon
would give `sum A_i^2+3 sum B_i^2=4` and `sum A_i*B_i=1`. This forces one
same-sign A_i=B_i=+1 or-1 at a single position, violating sum z_i=0. Thus
there is no hidden zero-height E8 enlargement.

## Direct interval and complete boxes

The proposed upper root has sum29+17s=(4+s)(5+3s), square2, and the two
nonzero old products -1 and-(1+s). Retaining n1 through n7 gives the stated
E8 finite diagram. Direct expansion yields

    Z0=11+8s,  T0=5+3s,
    Z0^2-8(s-1)T0^2=9,
    delta0=(4+10s)/71.

The six profile equations and all placements agree with finite extensions
of A7; the central placement4 is affine and correctly absent. Setting X=2T
gives D=2(s-1), with inverse(1+s)/4 exactly as implemented. From positive
X,Z and X/Z<=X0/Z0, the two direct squared bounds follow with inequality
in the code's direction. The conjugate equation is positive definite.

The enlarged ranges follow without a numerical cutoff: from
`0<Z<200, |Z'|<=12` one obtains
`-5<=A<=105,-3<=B<=70`; from `0<X<164, |X'|<6` one obtains
`-2<=C<=84,-1<=E<=56`. Boundary signs in `boxes()` and
`interval_solutions()` retain all possible roots. The exact square index
does not omit a field sign: positive X is required by the outgoing mirror
orientation, and all positive algebraic integers in the box are indexed.

## Candidates, first sector and actual endpoint

The reconstruction congruences, root norm and all eight old wall products
check. There is correctly **no** restriction p>=0 or z9>=0; the released
wall is tested as z9-p<=0. Because q is2 times a unit,2/q is integral, so
the M and ambient reflection matrices are integral. Nonprimitivity would
produce a unit-norm vector of M and is excluded by the mod2 argument.

For tied collisions, `T_L*r_R-T_R*r_L` has zero time, zero ninth coordinate
and zero first-eight sum. Its A7 coefficients are exactly the seven
nonnegative prefixes tested by `simple_dominates(L,R)`. Thus the selected
inequality implies all simultaneous candidate inequalities in the retained
sector; the mutual-domination assertion checks projective uniqueness.

The homogeneous endpoint belongs to M, is orthogonal to the retained A7
and selected root, and is strictly inside the released eighth wall. I
rederived its square as `-9(s-1)R<0`; its time is9Z>0. The source also
checks feasibility against every enumerated candidate. The first-endpoint
interpretation uses the written complete norm/profile inventory and the
locally finite reflection-chamber argument, not merely this finite Gram
calculation. Native angles with arbitrary unenumerated mirrors are not
separately censused by this worker.

The source binds both inputs and itself before/after, refuses overwrite,
bounds JSON, and replays through canonical JSON. Its45-second alarm does
not replace the root wrapper's shared lock, resource gates or memory cap;
it owns no inner global lock. Replay is the same implementation, not a
second independently implemented enumeration. No generation outcome is
asserted in this source review.

The proposed output concerns one opposite A7 edge of the full M chamber.
It proves neither compactness nor nonreflectivity of M or the original H9
host; those require an independent global argument or full cone control.
