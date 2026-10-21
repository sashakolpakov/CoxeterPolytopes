# TOOT: Exact marked bridge from the surface-free `T2` compiler to `C01`

Status: **datum-specific `IC0/R` bridge closed; existing `C/V` imported;
variable-class `U` open** (2026-08-08).

This certificate is separate from the Gram-only hostile precompiler.  The
precompiler remains an independent surface-free `IC0` calculation whose sole
mathematical input is the pinned 14-by-14 facet Gram.  Only after that archive
was frozen does this bridge verifier read the archived `(18,4,0)` lattice,
factory action and source walls, and existing `C01` strong-promotion archive.

It proves exact labelled recovery for this single `A18_4_0_C01` datum.  It is
not a new polytope, a class-wide `IC2/IC4` theorem, or forward universality for
a variable `H9/14` class.

## Exhaustive sheet-order search

The hostile fold-back has four fixed roots and ten exchanged root pairs.  The
factory has the same 14 active wall orbits and 24 short source walls.  The
verifier tries all `2^10=1024` orders of the two sheets.  Exactly four masks
make the ordered 24-root Gram equal to the hostile one:

```text
192, 440, 583, 831.
```

These differ only by compatible choices of sheet labels.  The frozen witness
uses mask `831`, with ordered factory source-wall indices

```text
1,9,3,8,4,7,5,6,17,11,16,12,15,13,
14,18,24,19,25,32,22,31,23,30.
```

## Literal integral marked isometry

Let `G_h,T_h,R_h` be the hostile rank-18 Gram, column-convention involution,
and ordered 18-by-24 root matrix.  Let `G_a,A_a,R_a` be the archived leading-
root Gram, factory row-convention action, and mask-831 ordered source-root
matrix.  The JSON records an explicit 18-by-18 matrix `P`, using

```text
x_archive = P*x_hostile.
```

It is integral, has determinant `+1`, and has integral inverse.  Direct exact
arithmetic verifies all three marked equalities:

```text
P^T G_a P = G_h,
P T_h       = A_a^T P,
P R_h       = R_a.
```

The frozen matrix digest is

```text
c1c7e49891a9467985e4df99d44ac41dff8bd21afa39c166dfc168a322ba1e14.
```

This is the required labelled integral recovery witness.  Agreement of
determinant, Smith form, discriminant invariants, action digest, or folded Gram
was not used as a substitute for these equalities.

The hostile primitive fixed and anti-invariant column lattices also transport
to the factory fixed lattice and the strong-promotion coinvariant lattice.
The two displayed basis changes are integral unimodular, of determinants
`-1` and `+1`, respectively.

## Full chamber and active-orbit checks

The audit does not stop at the folded Gram:

- the archived action is an exact order-two isometry of the rank-18 form;
- it permutes the complete set of 48 archived source walls, with an explicitly
  frozen involutive permutation;
- the full wall-norm histogram is `24*(-2)+24*(-4)`;
- all 24 source walls used by the active marking have square `-2`;
- the ordered source lists recover factory orbit indices
  `1,3,4,5,6,8,9,10,11,12,13,16,17,20` exactly;
- each singleton is fixed and each ordered pair is reversed by the action;
- each folded sum equals the factory
  `ambient_orbit_sum_primal_coordinates` row entry-for-entry; and
- the 14 transported folded sums reproduce the pinned facet Gram exactly.

Thus the isometry identifies the hostile fold-back with the actual archived
action and labelled active source-wall datum, not merely with an abstract
isometry class.

## Obligation boundary

The logical status is intentionally split by obligation.

`IC0` is closed for the single pinned facet Gram by the separate hostile
precompiler.  `R`, exact labelled integral recovery, is now closed for this
specific `C01` datum by `P` and the root/orbit equalities above.

Through that exact identification, the already-pinned strong-promotion result
supplies `C`: the factory cone is the complete invariant nef section, with all
restricted long walls redundant.  The pinned exact factory cone supplies
`V`: dimension nine, 14 irredundant facets, finite volume, noncompact, with
`90` finite and `20` ideal vertices.

`U` remains open: nothing here proves that every realization in a nontrivial
variable `H9/14` Coxeter class reaches this or another finite arithmetic
package and surface.  Consequently this datum bridge does not establish the
class-wide surface-softening equivalence.

## Reproduction

```bash
python3 scripts/audit_a18_c01_t2_marked_bridge_hostile.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_a18_c01_t2_marked_bridge_hostile_audit.py
```

Artifacts:

- [`a18_c01_t2_marked_bridge_hostile_audit.json`](../data/enumerations/a18_c01_t2_marked_bridge_hostile_audit.json)
- [`audit_a18_c01_t2_marked_bridge_hostile.py`](../scripts/audit_a18_c01_t2_marked_bridge_hostile.py)
- [`test_a18_c01_t2_marked_bridge_hostile_audit.py`](../tests/test_a18_c01_t2_marked_bridge_hostile_audit.py)

Frozen bridge JSON SHA-256:

`92396280c6b276979a005f9732408944c97fc71c08998d1f127833c7bacd6595`
