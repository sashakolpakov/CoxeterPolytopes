# Independent proof review: through-cusp walls cannot cut the H20 escape

**PASS:** every norm-two wall r acute to all 135 retained walls, with
Coxeter contacts and r.nu=0, satisfies r.E1<0. No hypothesis that r retain
the old interior point is used. Thus any finite collection of such walls
leaves this strict spacelike escape in place.

The theorem is local to the retained chamber and the specified null
direction. Walls with r.nu<0, deletion or replacement of old walls, and
other H20 constructions remain outside it. In particular it does not
rule out finite-volume completions using walls away from this cusp.

The [primary proof](fresh_h20_null_link_escape_gate_20261016.md) has two
independently checked layers:

- [Physical geometry reconstruction](../scripts/fresh_h20_null_link_escape_gate_20261016_geometry_review.py)
  and [coefficient certificate](fresh_h20_null_link_escape_gate_20261016_geometry_review.json).
- [Independent dual checker](../scripts/fresh_h20_null_link_escape_gate_20261016_review.py)
  and [complete case certificate](fresh_h20_null_link_escape_gate_20261016_review.json).

Generation and fresh saved replay pass for both layers. The reviewer
imports earlier independent field and interval engines, not producer
arithmetic or an optimizer. The optimizer's status is never evidence
for an infeasible or successful case.

## Physical and domain coverage

The geometry checker reconstructs the 107 old normalized roots from
literal algebraic coordinates, the 26 chain roots using the independent
recurrence, and both later caps from their exact recipes. It rebuilds
nu, E1 and the positive four-frame using the independently reviewed
projection construction. Norm and orthogonality equalities follow from
that construction's algebraic identities, not interval overlaps.

The [independent contact-domain theorem](fresh_h20_null_link_profiles_20261016_review.md)
covers 125 finite profiles and two continuous tails. For a finite profile,
the normal has the complete parameterization

```
r = p + sum z_i*f_i + lambda*nu,
p = -d^t G_W^-1 W,
||z||^2 = 2-p^2.
```

Here the 13 white roots W span A2+A5+A4+2A1. The positive affine G2
relation forces a through-nu acute wall to be right to that whole
triple. These conditions leave exactly four positive free coordinates
and one unrestricted null height. Replacing the norm sphere by its ball
is a necessary relaxation, including all 22 zero-radius profiles.

For a tail at root 17 or 31, the fifth bounded coordinate is d/sqrt(2),
the ball radius is sqrt(2), and sqrt(3/2)<=d/sqrt(2)<=sqrt(2). These are
whole continuous interval domains, including arbitrarily large finite
labels and the parallel endpoint.

The checker independently assembles the acute inequalities from all 119
old nonincident walls. At roots 2, 5, 61 and 87 it additionally takes
every one of the 16 choices between a right contact and a contact of
magnitude at least one. The right choice includes both inequality
directions. Other nonincident contacts remain merely acute. This covers
all 127*16=2,032 domains without imposing an arithmetic or label cutoff.

## Exact dual inequalities

For constraints A_z*z+A_h*lambda<=b and ||z||<=R, the objective is
e+c_z.z+c_h*lambda. Nonnegative multipliers satisfying the exact height
identity mu^t A_h=c_h give

```
objective <= e + mu^t b + R*||c_z-A_z^t mu||.
```

The independent checker takes the proposed rational multipliers as
untrusted witness data. It checks their indices and signs, reconstructs
all rows and bounds from its own physical coefficients, and defines one
pivot multiplier by the exact height equation. Rational outward bounds
prove that this multiplier is nonnegative. The remaining residual has
only bounded coordinates; its Euclidean norm receives an independent
rational upper bound. No small residual in the unbounded height is
accepted.

For infeasibility the target objective is zero. The same calculation
with mu^t A_h=0 gives a strict upper bound below zero on a quantity that
must be nonnegative in every feasible point, a contradiction.

| Case certificate | Count |
|---|---:|
|Exact infeasibility dual|1,739|
|Strict negative E1 objective bound|292|
|Separate homogeneous zero-line exclusion|1|

The weakest independent objective bound is approximately -0.157263820259;
the weakest infeasibility margin is approximately -0.000330346801.
The certificate stores rational bounds, and every bound is strictly
negative. These margins are not floating-point acceptance tolerances.

## The remaining homogeneous case

The zero profile and four right choices make r right to the full old
18-root support: its remaining member, root 1, is already right by the
affine G2 relation. Consequently r lies in the exact three-dimensional
Lorentz frame. The right conditions at roots 134 and 0 reduce it to a
single line. The reviewer checks that both denominators used to solve
these equations are nonzero, so the zero time-coordinate branch also
forces the entire vector to vanish.

The line generator has strictly negative pairing with old root 95 and
strictly positive pairing with old root 107. Acuteness to both forces
its scalar to be zero. This contradicts norm two and eliminates the
sole case where a homogeneous ball relaxation contains the origin.

The [previous independent cusp review](fresh_h20_affine_a17_completion_20261016_review.md)
proves E1 is a strict spacelike escape of the 140-wall chamber as well.
The theorem applies to any subsequent through-nu additions, since they
must in particular be acute to the retained 135 walls. A finite set of
such additions still has E1 strictly feasible and therefore cannot give
finite volume. A productive escape cutter must leave this cusp.

```sh
python3 scripts/fresh_h20_null_link_escape_gate_20261016_geometry_review.py --verify
python3 scripts/fresh_h20_null_link_escape_gate_20261016_review.py --verify
```
