# Exact congruence reduction from 51 to 34 H19 endpoint profiles

19 October 2026. New source-only mathematical packet. The preceding
`compact_h9_discriminant19_edge_pencil_obstructions_20261019.md`, SHA256
`41e78966ecb8137714347e79ddc5e434044409a85d5d98bcecd437350e8ace0a`,
is preserved unchanged. No scientific computation, finite census,
optimizer, Lean process or build ran.

**Results.** A complete modulo-3 argument leaves 40 of the preceding
52 Schur-positive labeled finite endpoint profiles (51 after its
all-zero exclusion). The coordinating agent's proposed prime-D
square condition is valid and removes six further profiles, leaving
34. Every retained scalar K is a square in Q(sqrt(5)); this is a
useful exact description, not an existence proof for their roots.
An additional observation of the finite57 agent also closes the
entire previously open omega3 short pencil modulo D.

None of these results excludes all short nonorthogonal endpoint roots
or proves whole-lattice nonreflectivity. The actual twenty-six-wall
seed is already known to have infinite volume by a separate escape;
a compatible completing mirror remains open.

## 1. Contact conventions and the complete scalar norm equation

Use phi²=phi+1, O=Z[phi], D=1+5phi, and the native form

    G=[[-2phi,phi],[phi,2]] + C,
    C=C_E8 with edges01,12,23,34,45,56,27.

The initial omega2-edge's eight white walls have positive-definite type

    W=A1(binary old0)+A2(E8 nodes0,1)
      +A4(E8 nodes3,4,5,6)+A1(E8 node7).

For a canonical short endpoint normal q=(x,y,z), let beta be the
nonnegative absolute product with old0. Let eta_j be its absolute
product with -e_j for j!=2. Each contact is 0,1,phi by the complete
finite-angle argument in the preceding report. Consequently

    phi*x+2y=beta,
    z=T*omega2 + sum_{j!=2} eta_j*omega_j,
    omega_j=C^{-1}e_j,
    R=sum_{j!=2} (C^{-1})_2j * eta_j,
    Q0=sum_{i,j!=2} eta_i*(C^{-1})_ij*eta_j.

Here T=(Cz)_2 is integral; all quantities are in O. There is still a
separate parity condition beta-phi*x in 2O for native y-integrality.
No congruence below silently discards that condition.

Since omega2²=30, the exact norm-two equation is

    beta²+2Q0+4TR+60T²-Dx²=4.                       (1)

Multiplying by15 and completing the square gives exactly

    A=30T+R,
    A²-15D*x²=K,
    K=60-15beta²-30Q0+R².                          (2)

This verifies the coordinating agent's identity with every factor
of2 and15 retained. Conversely (2) alone is insufficient: one also
needs A=R modulo30O, the parity condition, and all chamber inequalities.

## 2. A useful exact block identity for K

Writing f for the E8 indices other than2, the block inverse identity is

    (C^{-1})_ff
      -(C^{-1})_f2*(C^{-1})_2f/30 = (C_ff)^{-1}.

Thus 30Q0-R² is exactly 30 times the norm of the prescribed contacts
in the white A2+A4+A1 inverse Gram. Cross terms between its three
components vanish EXACTLY, not merely modulo a prime. Including the
binary A1 yields

    K=30*(2-eta_W^T*W^{-1}*eta_W).                 (3)

The previous Schur argument proves that at most one node of each
white component is contacted. Let a be the sole A2 contact, e the
sole A4 endpoint contact, i the sole A4 internal contact, and gamma
the E8 node7 contact; absent contacts have value0, and at most one
of e,i is nonzero. Equation (3) specializes to

    K=60-15beta²-20a²-24e²-36i²-15gamma².            (4)

The positions within the A2, and the two possible endpoint or
internal A4 positions, remain distinct labeled choices even though
they have the same scalar K. The two A1 are also distinct.

## 3. The full modulo-3 test: 52 profiles reduce to 40

The field O/(3)=F9 has phi²=phi+1. Its conjugation is phi ->1-phi,
and N(D)=-19=2 in F3. A nonzero square in F9 has norm1. Conversely
the four nonzero squares exhaust the four elements of norm1, so this
is an exact square test.

The complete omega2 row is

    (C^{-1})_2j=(10,20,30,24,18,12,6,15).

If the A2 contact is nonzero, R=a or2a modulo3 according to its
chosen node. Hence R is nonzero in F9. In (1) the free T coefficient
4R is then nonzero, so for any chosen x there is a unique T satisfying
the reduced norm equation. This particular modulo-3 test excludes
none of the A2-contact profiles. It makes no claim about other primes,
parity, or geometric feasibility.

Suppose instead that neither A2 node is contacted. Then R=0 modulo3.
The required inverse entries satisfy

    (C^{-1})_33=(C^{-1})_66=2 mod3,
    (C^{-1})_44=(C^{-1})_55=0 mod3,
    (C^{-1})_77=2 mod3.

All A4/node7 cross entries are0 modulo3; there is at most one A4
contact. Reducing (1) therefore gives exactly

    D*x²=beta²+gamma²+e²-1 in F9.                 (5)

An internal A4 contact contributes nothing to (5), but its presence
has not been discarded from the Schur-positive profile.

For profiles with only contacts1, put k=beta²+gamma²+e², the integer
number of contacted nodes among these three specified possibilities.
Here k is0,1,2 or3. If k=1, equation(5) permits x=0 modulo3. For each
other k, the nonzero right side divided by D has norm2, so is not a
square. Thus EXACTLY ONE of those three contacts must be present.

Of the 18 Schur-positive zero/1 profiles without A2, precisely eight
survive: A4 endpoint alone (2), one A1 alone (2), and A4 internal
plus one A1 (4). The other ten, including the all-zero profile, fail.

For profiles involving a contact phi but no A2, the prior Schur
list contains only a single A1 contact phi, optionally with contact1
at the other A1. In the first case the right side of (5) divided by
D is phi/D=phi^{-2} modulo3, a square. In the second it is phi²/D,
whose norm is2, so it is not a square. This removes two further
labeled profiles. The remaining profiles involving A2 were already
covered by the nonzero-R argument.

The exact accounting is therefore 52-10-2=40. Compared with the
preceding report's 51 nonzero profiles, this is eleven new exclusions.

## 4. The prime-D condition: 40 profiles reduce to 34

The ideal(D) has norm19; its residue field is F19 with phi=15.
Reducing (2) gives the necessary condition

    A²=K mod(D).

Thus K must be a square modulo19. The squares, including0, are

    {0,1,4,5,6,7,9,11,16,17}.

From (4), since phi²=16 modulo19,

    K=3+4beta²-a²+14e²+2i²+4gamma² mod19.           (6)

Here is the COMPLETE list of40 profiles after Section3. A label
"3" means contact1; "5" means contactphi, matching the corresponding
Coxeter angle. Multiplicities retain the actual white-node labels.

|Contacted components|Multiplicity|Exact K|K modulo19|Passes prime-D test|
|---|---:|---:|---:|---|
|A2(3) alone|2|40|2|No|
|A2(3)+A4 endpoint(3)|4|16|16|Yes|
|A2(3)+A4 internal(3)|4|4|4|Yes|
|A2(3)+one A1(3)|4|25|6|Yes|
|A2(3)+A4 endpoint(3)+one A1(3)|8|1|1|Yes|
|A2(3)+both A1(3)|2|10|10|No|
|A4 endpoint(3) alone|2|36|17|Yes|
|One A1(3) alone|2|45|7|Yes|
|A4 internal(3)+one A1(3)|4|9|9|Yes|
|A2(5) alone|2|40-20phi|6|Yes|
|One A1(5) alone|2|45-15phi|10|No|
|One A1(5)+A2(3)|4|25-15phi|9|Yes|

Exactly six profiles fail, giving 34 survivors. This is a necessary
local filter, not a proof that the retained congruences, integral
norm equations, or actual chamber inequalities are simultaneously
solvable.

## 5. Every retained K is a field square

The retained rational values are

    1,4,9,16,25,36,45,

all squares in K, including 45=(3sqrt(5))². The two remaining values
also have explicit square roots:

    40-20phi =20phi^{-2} =(2sqrt(5)*phi^{-1})²,
    25-15phi =5phi^{-4} =(sqrt(5)*phi^{-2})².

Thus every remaining branch has a generalized Pell equation with
square right side. This may support a relative-unit or residue-orbit
approach. It does not prove existence of a nonzero-x root satisfying
A=R modulo30O, beta-phi*x in2O, T>=0, all other acuteness and
angle requirements, and an actual closing position. Those obligations
are not replaced by the trivial field solution x=0, A=sqrt(K).

## 6. Additional all-height obstruction for the omega3 short pencil

The finite57 agent observed the following after its bounded pencil
test; I independently checked the argument on paper. A canonical
short root in that orthogonal pencil would satisfy

    10t²-Du²=1.

Modulo(D), this forces t²=2 in F19, a nonsquare. Consequently this
entire pencil has no canonical short root at any height. The preceding
report already excludes its long roots modulo sqrt(5), so both norm
classes in the pure omega3 pencil are now obstructed. Its bounded
empty search is not the basis of this all-height conclusion.

The omega2 edge can still be closed, if at all in an append-only
Coxeter construction, by a short nonorthogonal root in one of the
34 necessary labeled profiles. No such root is constructed or
excluded here. No whole-lattice nonreflectivity or new finite-volume
example follows from these local arithmetic reductions.
