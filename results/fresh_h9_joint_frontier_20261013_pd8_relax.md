# Continuous relaxation for old-edge endpoint profiles

This is an analytical formulation and a **numerical dual-proposal diagnostic**, not a completed exclusion theorem. It concerns caps meeting one of the 302 old open edge supports of the fixed 27-wall chamber: 297 supports with a finite endpoint and five vertex-free lines. A cap meeting a rank-seven face instead can evade this mechanism; the constructive search is now pursuing such a cap.

For a spherical endpoint profile, normalize every wall and the proposed cap to norm two. The eight prescribed pairings determine a projection `P` onto the positive-definite support space. Its residual norm is `rho=2-P²>0`. Choose an orthonormal basis `u,e` of the orthogonal Lorentz plane, with `u²=-1`, `e²=1`. Every residual root has one of the two forms

```text
r(t) = P + sqrt(rho) [sinh(t) u + sigma cosh(t) e],  sigma=+1 or -1.
```

Set `T=exp(t)>0`. Multiplying each inequality by positive `T` turns every old-wall acuteness condition `(r,w)<=0` and the escape-cut condition `(r,E)>0` into a quadratic polynomial inequality in `T`. A complete continuous test therefore intersects finitely many algebraic intervals on both branches; it must retain all closed acuteness endpoints and strict escape-cut boundaries. No finite sampling of ninth anchor labels establishes this completeness. A rigorous implementation can isolate the quadratic roots and compare them, or certify interval covers with exact outward coefficient bounds.

An isolated `I2(m)` profile has one contact `-c`, with `c=2cos(pi/m)` and `1<=c<2`, projection `P=-c*w/2`, and residual norm `2-c²/2`. Allowing the whole closed interval `1<=c<=2` gives a continuous superset of every integer label. Excluding that superset would settle the entire tail. A surviving continuous interval would require testing the discrete cosine sequence; it must not be reported as a reflective cap without the remaining angle gates.

Before the norm equation, a simpler certificate may suffice. Find

```text
E = lambda1*w1 + lambda2*w2 + sum(mu_i*s_i),  lambda1,lambda2>=0,
```

where the `s_i` are the eight support normals. A profile with `(r,s_i)=-h_i` then satisfies `(r,E)<=-sum(mu_i*h_i)`. If the last sum is nonnegative, the cap cannot cut `E`, without any bound on its hyperbola parameter. The [proposal script](../scripts/fresh_h9_joint_frontier_20261013_pd8_relax.py) searches these two-wall identities using floating-point linear algebra only. It uses the exact integer escape recorded by the [old-bank gate](fresh_h9_makarov_joint_completion_20261013_old_bank_gate.json), not the later numerical relaxed-bank gap.

Across 10,200 finite profile records and 432 symbolic tails, the [numerical output](fresh_h9_joint_frontier_20261013_pd8_relax.json) gives:

| Diagnostic | Records |
| --- | ---: |
| Proposed strict dual bound | 7,732 |
| Proposed zero dual bound | 253 |
| A dual exists but the continuous norm gate is still needed | 1,203 |
| No proposed two-wall dual | 1,444 |

These are unverified proposals. No exact exclusion is asserted for any of these counts, and absence of a two-wall certificate does not prove a cutter exists. The newly found rank-seven construction takes priority over exactifying this restricted old-edge test.

For that constructive continuation, the [exact PD7 handoff](fresh_h9_joint_frontier_20261013_pd7.json) exhausts 2,234 finite-pair seven-cliques and finds 1,339 positive-definite two-face supports. They have 1, 2, 3, 4, or 5 existing edges with respective frequencies 20, 671, 483, 68, and 97. Its [checker](../scripts/fresh_h9_joint_frontier_20261013_pd7.py) includes every support and its certified edge incidences. This is geometry of the current chamber, not an assertion that any face has a compatible completion.
