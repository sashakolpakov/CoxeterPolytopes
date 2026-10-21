# Review of universal small-A component energy

20 October 2026. Source and paper review only; no producer, checker,
optimizer, Lean or TeX process was run.

Reviewed in full:
`results/compact_universal_small_a_component_energy_20261020.md`,
SHA256 `801581c4b855289aca7982b2d0f4bb6f5849f151a78f40b86c9d403e84f415e6`.

**PASS within the stated scope.** At an actual vertex farthest from a
strict interior point of a compact hyperbolic Coxeter polytope, an A3
component costs more than 1/4 and an A4 component more than 2/7, with
arbitrary finite spectator components. A whole vertex forest consisting
of A1, A2, A3 and A4 factors is excluded in every dimension at least 14.
This is not a dimension bound for unrestricted vertex forests.

I checked the following substantive interfaces.

* A finite retained spectator has at most one contact with an outside
  wall. Its positive projection cost is at least 1/4, and at least 1/3
  for a nontrivial factor. Together with the actual cross-component
  opposite obstruction this proves small-height isolation without a
  spectator catalogue.
* The pair/logarithmic-concavity argument forces one actual opposite
  wall under the proposed A3 and A4 energy thresholds. The native
  multi-contact exclusions and the four single-contact inequalities
  then yield the stated universal component bounds. The A4 asymmetric
  endpoint minimum permits either endpoint to carry the larger bond.
* The short A1 classification exhausts zero, one and multiple spectator
  possibilities. Only isolated dotted I and dotted-plus-A1_3 J remain;
  the interior inequality gives J height square greater than 1/3.
  Actual-wall orthogonality then permits at most one short A1 direction,
  and the light J coordinate funds its heavy-coordinate deficit.
* The short A2 inventory includes the finite isolated labels at least
  7, isolated dotted contacts, and exactly the named A2_m+A1_3 cases.
  Actual mate recognition and distinct-wall shadows are used in their
  proper branches. Only isolated 7 can debit the A2 lower bound.
* Both propagation bridges hold for actual opposites in other whole
  components. In particular their projection-cost bounds use q<1 from
  a retained finite white component, and their directional derivative
  bound does not omit the other spectator costs. The isolated-7 and
  small-I cases force enough remaining direction energy to exceed one.
* The final coefficients and the n=14 parity argument in (11) are exact:
  when the isolated I wall is present, a=1 forces an odd positive number
  of A3 factors, whose surplus exceeds the remaining A1 debit. Higher
  dimensions and a>=2 give the stated contradiction directly.

The inputs remain the written actual-farthest-point, reciprocal-shadow
and cross-component-opposite lemmas, all explicitly identified in the
source. I previously co-developed these mechanisms with the author and
root; this is an outside-author full-source audit, not an independent
discovery or a claim of formal verification. No new flaw was found.
