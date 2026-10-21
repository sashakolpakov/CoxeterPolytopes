# Correction: all three long profiles on the new H19 edge are impossible

19 October 2026. Independent written coordinate and arithmetic review only.
No checker, matrix calculation program, search, enumeration, Lean or build
was run. The reviewed report and every pinned source remain unchanged.

**Verdict.** The16-profile classification in
`h19_next_edge_sixteen_profiles_research_20261019.md` reduces to **11 short
profiles and no long profiles**, not12 profiles. Its long-root equations(3)
and(4) are correct, but their reduction modulo the discriminant prime is
incorrect. The reported relation phi^5=1 modulo19 is also incorrect:
phi^5=2. This removes the proposed surviving long contact-b branch at every
height; it does not settle any of the remaining11 short profiles or the
possibility of completing the chamber.

## Re-derivation from the native form and reflection integrality

Let O=Z[phi], phi²=phi+1, b=4+phi and D=phi*b=1+5phi. The native form in
the26-wall checker is

    G = [[-2phi,phi],[phi,2]] + C_E8.

The old0 white root is(0,−1;0), and old3 through old8 are the six negative
E8 simple coordinate roots with indices2 through7. Thus the seven short
white orthogonality conditions on a candidate q=(x,y;v) are

    phi*x+2y=0,       (C_E8*v)_j=0 for j=2,...,7.

The inherited complete primitive norm alphabet is{2,2b}, up to unit
rescaling. For a canonical integral long reflector q²=2b, the integral
reflection criterion is Gq/b in O^10. In particular C_E8*v belongs to
bO^8. Since phi is a unit, phi*x+2y=0 with x,y integral implies x=2u and
y=−phi*u for u in O. The two free E8 covector entries can be written
b*tau0,b*tau1 with tau0,tau1 in O. Unimodularity gives exactly

    q=(2u,−phi*u; b*(tau0*omega0+tau1*omega1)),
    omega_j=C_E8^(-1)e_j.

The actual checker displays

    omega0=(4,7,10,8,6,4,2,5),
    omega1=(7,14,20,16,12,8,4,10).

Hence omega0²=4, omega0·omega1=7 and omega1²=14. The binary contribution
to q² is−2D*u². Dividing the full identity q²=2b by2b yields

    −phi*u²+b*(2tau0²+7tau0*tau1+7tau1²)=1.       (A)

Thus equation(3) follows directly from the native coordinates; no
congruence of the reported completed equation is being assumed.

## The actual z wall fixes the second equation

The26-wall checker defines

    z=(2phi^4,−phi^5; b*phi²*omega0),
    Gz=(−b*phi^5,0; b*phi²*e0),       z²=2b.

For a finite endpoint, the inherited native finite-angle alphabet between
two canonical long roots gives q·z=−ab with a=0,1 or phi. This is exactly
the right/angle3/angle5 list in the16-profile report. Since
v0=b*(4tau0+7tau1), direct pairing with this **actual** covector gives

    q·z=b*phi²*(b*(4tau0+7tau1)−2phi³*u)=−ab.

Dividing by the nonzero b*phi² proves

    b*(4tau0+7tau1)=2phi³*u−a*phi^(-2).          (B)

This verifies equation(4), including both powers of phi and the sign of
the prescribed acute contact. It does not infer an unspecified wall
contact from the report's abstract profile alone.

## Correct reduction in O/(b)

The ideal(b) has norm19; the quotient is F19 with phi=−4=15. Since phi is
a unit, (b)=(D). The needed residues are

    phi²=16, phi³=12, phi^4=9, phi^5=2,
    2phi³=5, phi^(-2)=6.

Equation(B) therefore becomes5u−6a=0, so

    u=5a modulo(b),

because5^(-1)=4 and6*4=5 modulo19. Equation(A) independently reduces to
−15u²=1. Substituting u=5a gives

    5a²=1,  equivalently a²=4 modulo19.

The three allowed contact parameters a=0,1,phi have residues0,1,15 and
squares0,1,16. All three contradict a²=4. This uses a necessary condition
on every canonical integral long reflector, so it excludes all three
profiles at every height without a bounded search or further wall tests.

The report's statement u=a/2 would be u=10a in this quotient and is the
incorrect step. Its subsequent surviving a=1 branch and proposed long
search domain must consequently not be used.

## Remaining scope and exact source bindings

The short-root count13 is unchanged. Its two no-A6-contact profiles with
beta=0 or1 are excluded by the report's original modulo7 norm equation:
7 is inert, N(phi^(-3))=−1, and the forced square has nonsquare norm−1 in
F7. Those exclusions leave11 short profiles. This review establishes no
new short-root exclusion and no compatible mirror, finite-volume chamber,
whole-lattice nonreflectivity or dimension result.

The complete norm alphabet and finite-angle alphabet remain explicit
inherited arithmetic inputs. Only the long branch's coordinate reduction
and congruence correction were independently rederived here.

* Reviewed16-profile report: SHA256
  `6dc87d20c5cf207d90aabe28d3f1aa1b8c38596484f8fd6ab98c6532d7ff95cb`.
* `scripts/check_h19_twenty_six_wall_seed_20261019.py`: SHA256
  `f0ac1d28225a37be17e9b61f2e436ea2fc3b1511efb32b02212b33539fecb363`.
  Its native form, omega vectors, z coordinates/covector and norm2b
  reflection-integrality condition were read; `check()` was not invoked.
* Native prefix `data/candidates/h19_vinberg_prefix.json`: pinned SHA256
  `91c16b581e5763d6a283b6818f0e087e941fbb987b6b05a7e283febf1e32d0f1`.
  Its literal height-zero simple roots fix the white-coordinate mapping.
* Complete norm-alphabet source
  `compact_h9_discriminant19_mixed_norm_seed_20261019.md`: inherited
  SHA256 `70fa010732d837d5c044701b4c18ed5e68cb36c18cc0954dc7d7f369250a5c87`.

The earlier16->12 report is preserved as historical evidence. This
separate correction establishes16->11 for the same actual edge and
canonical integral norm classes, conditional on the stated inherited
profile/angle completeness inputs.
