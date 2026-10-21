# Independent review: full-reflection volume bridge

20 October 2026. Paper/source review only; no scientific execution.

**PASS within the stated conditional structural scope.** I read the full
addendum
`results/bugaenko_facet_full_reflection_volume_bridge_addendum_20261020.md`,
SHA256 `e61a5531051aac247382b9e1ec2f4536bb7a9b89fba3352febf7fe1bf01cf6a3`.
It fills the preceding volume-upper-bound interface for full integral
reflection chambers, with ineffective constants. It supplies no new
dimension exclusion or additional host exclusion.

For the full integral group `Gamma=O^+(R)`, conjugation preserves the set
of all integral reflections. Thus their subgroup `W` is normal, and its
finite covolume gives finite index. The action on actual chambers gives
`Gamma=W semidirect A`, with `A` the finite chamber stabilizer. This
argument uses the full integral group, not maximality of `W` among every
reflection group in the ambient real isometry group. Admissibility over
the nonrational golden field makes the underlying form anisotropic, so
the finite-index chamber is compact.

I checked the elementary reduction-modulo-3 bound separately. Since 3
is inert, `O/(3)=F9`. In a prime-order kernel element `I+3^j B`, the
first binomial term has strictly smaller 3-adic valuation than every
later term, including order three. Hence the congruence kernel is
torsion-free on finite subgroups, and `|A|<3^200` follows. The full
integral group itself is congruence; the proof does not require `W`
to be congruence.

I directly checked the analytic inputs in the primary
[Agol--Belolipetsky--Storm--Whyte paper](https://homepages.math.uic.edu/~kwhyte/ahfinite.pdf):
Theorem 3.1 and its proof, printed pp. 4--5; Theorem 4.2, p. 6;
Theorem 5.9, p. 11; and the quotient application in Section 6,
pp. 12--13. The underlying piecewise conformal framework allows the
mirrored orbifold boundary. The published application uses the same
finite chamber stabilizer and congruence ambient quotient interface;
their maximality step is replaced here by the preceding direct proofs.

The displayed constants check: fewer than `50m^2` cutting hyperplanes
give at most `20*50^9*m^18` pieces. Applying the spectral gap `15/4`
in source dimension nine contributes `(12/5)^(9/2)`, and passing from
the ambient quotient to the chamber contributes one further factor
`m`. The resulting exponent is therefore 19 and the stated uniform
upper bound with `3^3800` is valid.

The inherited collar lower bound has coefficient
`kappa=phi^4*pi^4/(24*|W(E8)|)`. Substituting
`q=4*cosh(ell)^2` in the earlier determinant bounds gives respectively
`16*sinh(2ell)^2` and `256*sinh(2ell)^2`, exactly as claimed.
The odd-glue basis change is integral because `phi^3=1 mod 2`; it keeps
the cross pairing one and multiplies the determinant by `phi^(6j)`.
The direct-glue change multiplies it by `phi^(2j)`. Consequently a
bounded absolute determinant norm gives bounded embeddings after the
stated reduction, hence finitely many integral coefficient pairs.

The preceding collar source was independently reviewed at
`results/review_bugaenko_facet_binary_collar_batch_guard_20261020.md`,
SHA256 `81f7aa5744b1c0cb1f00abb5542afc9f25cc1463236dc0831d7858d3d28d3e0e`.
The addendum correctly keeps full-chamber volume separate from the
volume of an arbitrary finite-index reflection subgroup. Its huge
finite determinant bound is not a feasible census and is not a
uniform dimension improvement. No source changes are requested.
