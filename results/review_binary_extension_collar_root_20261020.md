# Root review: binary reduction and uniform E8 collar

20 October 2026. Full written review of
`bugaenko_facet_binary_extension_and_collar_audit_20261020.md`,
SHA256 `774e87f7abb33680a06879a17a5be9f06e2218bf12a95023bcff5c7046ba1227`.
PASS as a structural research argument, not a dimension exclusion.
No scientific execution or new formal verification was needed.

The connected 3/5 facet graph permits normalization over the ambient
field. Any extra embedding fixing phi would retain an indefinite
rank-nine subform, proving the asserted exact field. Finite-index
commensurability makes each group element's eigenvalues integral; the
trace formula for two norm-two reflections then gives integral orbit
pairings. A finite coset decomposition bounds denominators of the whole
orbit. Cocompactness forces rank ten. Unimodularity splits E8, evenness
makes the time vector primitive, and reduction modulo2 gives exactly the
two displayed glue forms. None of these steps presupposes that the
original lattice was even or maximal.

For a primitive root of norm2d, reflection integrality and E8
unimodularity give m_E8=d*z. The quantity T=d*z^2/2 is integral and
has conjugate between0 and1. Its nonzero integer norm forces T>=1;
T=1 implies a pure E8 root. In the mixed case write T=A+B*phi:
B>=1 and A>=1, so T>=phi^2. Reverse Lorentz Cauchy proves the claimed
uniform distance, and the Fermi-coordinate integral gives phi^4/8.

One implicit classification detail also checks directly. A primitive
pure E8 root has unit d; a unit rescaling gives square2. Write that
vector a+phi*b with a,b in the ordinary even E8 lattice. Square2 forces
a^2+b^2=2 and 2(a,b)+b^2=0. Evenness and minimum norm2 imply b=0 and
a is an ordinary E8 root. Thus the angular sector used in the collar
really is the ordinary E8 Weyl sector.

In the odd glue case q and q-4 have the stated binary Gram expressions.
The square factor in q(q-4)=D*s^2 is integral and nonzero; the conjugate
bound q'(4-q')<=4 proves |N(D)|<=4q(q-4). The direct case gives the
weaker but valid constant64 for the determinant norm. These concern
the actual next binary wall, with no chosen norm-class assumption.

The infinite module family, explicit basis changes, seed escape and
limited unit-denominator point-inversion mechanism also check. The note
correctly distinguishes an integral lattice involution from a chamber
stabilizer. No uniform chamber-volume upper bound, exhaustive binary
classification, compact H9 example or improved ceiling follows yet.
