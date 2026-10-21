# D18 centre scalar: prepared logic and assumptions

`formal/CoxeterProofs/D18CentreScalarKernel.lean` is **prepared, uncompiled
source**. It requests five axiom prints and has **zero actual outputs**.
No worker, import, test, compiler, Lean or TeX job ran during preparation.
The module imports only `Std` and changes no existing successful kernel.
Root owns any later bounded compilation and separate actual axiom audit.

The mathematical application is Section5 of
[the complete D18+actual-I supporting proof](../results/compact19_d18_isolated_i_complete_centre_and_head_exclusion_candidate_20261021.md),
source SHA-256
`8b8e4b90f4059a85778fe0a9c033d80c10d02a932d877b0d993dc588973bed73`.
This module is intended to check only its scalar algebra, conditionally on the
application bridges below. It is not a formalization of that whole
source, a D18 exclusion or a dimension theorem. The source's root
mathematical review is distinct from this kernel's pending source and
execution reviews.

## Derived statements

The type R is an arbitrary linearly ordered commutative ring with the
same `Lean.Grind` classes used by the existing Std kernels. It can be
instantiated with genuine real geometry; no integer or rational proxy
replaces R. Every denominator below is cleared in the Lean statements.

1. `quadratic_completion` proves the exact identity for
   `q(z)=96874800*z²-1528000*z+6100`:

   ```text
   387499200*q(z)=(193749600*z-1528000)²+28961120000.
   ```

2. `quadratic_positive` obtains square nonnegativity from the ordered-ring
   laws and derives `q(z)>0` for every z. Neither this positivity nor
   the square's sign is an application premise.
3. `continuous_scalar_lower_bound` takes only the analytic lower bound
   `1000000*F>=34196500-1528000*z+96874800*z²`, and derives
   `10000*F>341904`. It does not assume the uniform scalar conclusion.
4. `k10_margin` derives `T>1` from that analytic bound,
   `100000*base>=87866`, `1000*f<=2731`, and
   `49*T>=49*base+F-14*f+10`. The intermediate exact inequality is
   `4900000*T>4901074`.
5. `k11_margin` uses instead `100000*base>=90246`,
   `10000*f<=28854`, and the same payment with final constant11.
   Its intermediate exact inequality is `4900000*T>4901534`.

Thus the respective cleared lower bounds on49T are49.01074 and49.01534.
The positivity and strict final margins are consequences, not supplied
dimension exclusions or norm-budget contradictions.

## Explicit application inputs

For the actual no-own finite-fork cap of Section5, the written geometry
defines a real parameter `.81<b<1` and

```text
x=1/b, z=x-7/6,
G=sqrt(18-2/b)-sqrt(-8+14/b-6/b²),
R=26-16/b+6/b², F=14*G-R.
```

The source derives an actual feasible future negative-core point using
positive combinations of actual wall normals. It checks all original
wall inequalities and computes the norm, time, displacement and pressure.
It then proves the radical chord/tangent bounds and the rational
comparisons `2.4494<sqrt6<2.4495`. Those analytic and geometric steps
produce the Lean hypothesis `analytic_lower`; they are **not** proved
by this module. No interval, square root or division operation appears
in its ring statements.

The source also supplies the complete prior local inventory, a genuinely
strong K10 or K11, the original farthest vector, positive reference
coefficients and all actual-foot pressures. Paying the supported endpoint
coefficient before adding the core pressure yields `paid_reference` with
T equal to the squared norm of that SAME original vector. This bridge
includes the full norm saving, the mixed support terms and the legal
initial/I adjustment. It is not an arbitrary assertion that T>1.

For strong5 the written branch has `base=qK+joint_gain`, where
`q10>.8244`, `q11>.8482`, and the once-counted initial-five/I gain
is>.05426. This gives the two base floors. For weak5 the stronger
initial references L10>.8870 and L11>.9104 also satisfy those floors,
without reusing I energy. The separate bounds f10<2.731 and
f11<2.8854 supply `height_upper`. None of these reference, radical
or geometric bridge claims earns Lean credit merely by appearing here.

Finally compact farthestness gives T<1. The ordered-ring conclusions
T>1 therefore contradict the application, but the kernel does not
construct the polytope model or prove its norm budget, exhaustive cap
inventory, other centre/head cases, arithmetic facts, or published
geometry. No global bound, new example, runtime success or axiom-audit
result is claimed while this source remains uncompiled.
