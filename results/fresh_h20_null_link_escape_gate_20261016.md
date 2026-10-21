# A complete through-cusp escape gate for the retained H20 chamber

Every further acute Coxeter wall through the specific future null ray
nu misses the strict escape E1 of the retained 135-wall chamber. Exact
generation has checked all 2,032 domains: 1,739 infeasibility
certificates, 292 strict escape upper bounds and one direct zero-line
exclusion. There are no unresolved domains. Generation and a fresh
optimizer-free saved replay both pass. Independent review of this new
gate is pending.

The conclusion is deliberately local: any norm-2 wall r acute
to all 135 retained walls, with `r.nu=0` and Coxeter contacts, satisfies
`r.E1<0`. It does not exclude walls with `r.nu<0`, changed old-wall sets,
or other H20 constructions. No condition that r preserve the old center
is imposed in this proof.

## Exact domain and coordinates

The [135-wall certificate](fresh_h20_complementary_cap_20261016.md) defines
the strict spacelike escape E1 and the affine G2 triple `(134,0,1)` with
positive null relation

```
nu = sqrt(3)*n134 + 2*n0 + n1.
```

If an acute wall r is through nu, every summand in `r.nu` is nonpositive,
so all three contacts with that triple must be right. The other old
incident walls form the positive block

```
A2[3,4] + A5[6,7,8,9,10] + A4[12,13,14,15] + A1[17] + A1[31].
```

The [complete contact-domain certificate](fresh_h20_null_link_profiles_20261016.md)
and its [independent review](fresh_h20_null_link_profiles_20261016_review.md)
give 125 explicit finite profiles and two continuous single-root tails.
The latter allow magnitude d in `[sqrt(3),2]` at root17 or root31, with
all other incident contacts right. They cover every remaining integer
label and the parallel endpoint. There is no all-label sampling cutoff.

Use the exact orthonormal four-frame f_i from the
[A17 cusp construction](fresh_h20_affine_a17_completion_20261016.md),
which has a separate [independent review](fresh_h20_affine_a17_completion_20261016_review.md).
For a fixed finite contact profile d on the 13 white roots W, put
`p=-d^t G_W^-1 W`. Every candidate has a unique expression

```
r = p + sum_(i=1)^4 z_i*f_i + lambda*nu,
r^2 = p^2 + ||z||^2.
```

The exact identities hold because f_i are orthonormal and orthogonal to
W and nu, and because the G2-orthogonal Euclidean quotient has dimension
17=13+4. Lambda is unrestricted. Replacing norm equality by
`||z||<=sqrt(2-p^2)` only enlarges the candidate set.

For either continuous tail use

```
r = sum z_i*f_i - (d/2)*n_j + lambda*nu,
||(z_1,...,z_4,d/sqrt(2))|| <= sqrt(2),
sqrt(3/2) <= d/sqrt(2) <= sqrt(2).
```

These are continuous whole-interval domains, not endpoint tests.

## The linear contact relaxation

All 119 old nonincident walls impose their acute inequality. Four of
them, roots2,5,61,87, are the remaining nonincident members of the old
positive 18-root support. At each of these four, a Coxeter contact is
either right or has magnitude at least1. The checker takes all 16
choices. A nonright choice uses only `r.n_j<=-1`; a right choice adds
the reverse of the already imposed acute inequality.

Thus the finite profiles have 123 linear inequalities in four bounded
coordinates and one unbounded height. The tails have two further
interval inequalities and one further bounded coordinate. The
optimization considers every one of the 127 times16 domains. All other
nonincident contacts are merely acute, so the relaxation includes every
actual Coxeter possibility, even with arbitrary dotted distances.

The [geometry source](../scripts/fresh_h20_null_link_escape_gate_20261016_geometry.py)
and [coefficient data](fresh_h20_null_link_escape_gate_20261016_geometry.json)
reconstruct the literal old normals and evaluate rigorous rational
enclosures for each linear coefficient. Exact Gram identities are
inherited from the explicitly defined frame, not inferred from intervals
that happen to contain zero. Generation and fresh saved replay of the
coefficient package pass.

## What a dual certificate proves

Write a domain as

```
A_z*z + A_h*lambda <= b,   ||z||<=R,
r.E1 = e + c_z.z + c_h*lambda.
```

For nonnegative multipliers mu with the exact height balance
`mu^t A_h=c_h`, every feasible vector satisfies

```
r.E1 <= e + mu^t b + R*||c_z-A_z^t mu||.
```

A strict negative upper bound therefore proves that the domain misses
E1. Alternatively, nonnegative multipliers with `mu^t A_h=0` and

```
mu^t b + R*||A_z^t mu|| < 0
```

prove that the domain is empty.

The optimizer only proposes multipliers. All but one are rounded to
nonnegative rational numbers. A recorded pivot multiplier is **defined
exactly** by solving the height balance; its nonnegativity and the strict
bound are then verified with outward rational interval arithmetic.
This exact defining equation removes any unbounded-height error from
floating-point stationarity. Positive rational row scaling does not
change the inequalities.

The [main checker](../scripts/fresh_h20_null_link_escape_gate_20261016.py)
replays these certificates without importing an optimizer. Its
[result](fresh_h20_null_link_escape_gate_20261016.json) records every
domain, all rational multipliers, the balance pivot and a strict bound.
Any unresolved domain remains explicit and prevents the final PASS
status; optimizer infeasibility alone is never accepted as proof.

## The sole homogeneous zero-profile case

When all 13 incident contacts and all four selected nonincident contacts
are right, r is right to the entire old 18-root support. It therefore
lies in the three-dimensional Lorentz frame. Right contacts with
wall134 and wall0 restrict it to one line. A generator with time
coordinate1 is defined exactly by

```
z = y134_time/y134_f,
e = -(P0_time + P0_f*z)/P0_e.
```

The coefficient checker proves its pairing with old95 is negative and
with old107 positive. Acuteness to both forces its scalar coefficient
to be zero. This excludes every norm-2 wall in this domain directly,
without a strict numerical margin at the zero vector.

The weakest certified escape upper bound is below -0.15726382024. The
weakest infeasibility inequality is below -0.00033034680 (after the
recorded positive row scaling). Both are strict rational interval
bounds, not rounded optimizer values.

Thus E1 stays strict after any finite coherent set of through-nu walls.
Such a set cannot produce a
finite-volume completion. In particular, completing the local cusp
link—although achieved exactly by the five A17 connectors—does not
remove this distinct open region. A further E1 cutter must have
`r.nu<0` while retaining the old affine G2 triple.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_null_link_escape_gate_20261016_geometry.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_null_link_escape_gate_20261016.py --verify
```
