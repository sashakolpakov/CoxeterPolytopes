# Prime811: complete next-edge domain and a firstness certificate contract

21 October 2026. Independent manual source review and one new row-reduction
lemma. No root enumeration, scientific import, checker execution, Lean or
TeX. No new polytope or literal-host exclusion is claimed. This is a
supporting construction audit, not a seventh manuscript.

## Reviewed source and verdict

Reviewed `compact_h9_prime811_construction_route_and_edge_certificate_20261021.md`,
SHA256 `f814a811e96d0519e2cb1c812370c6af73fcf705f0631d115407b8e4babef825`
(8812 bytes), under `results/`. Its six displayed inherited source pins
were independently streamed and all match. Their mathematical theorems
remain inherited inputs; this review rederives the new k5 arithmetic,
edge profiles, bounds and conditional isometry rather than reopening the
full Bugaenko calibration. I did not author the reviewed candidate.

The candidate is mathematically consistent at its stated CONDITIONAL
scope. In particular, it supplies a finite next-edge problem, not the
firstness conclusion. The best next step is one small exact certificate
for this entire bounded edge. There is presently no demonstrated earlier
wall, nor a proof that its displayed upper is first.

## Checks of the new arithmetic and completeness

With phi^2=phi+1, the data are

    u=3+5phi, a=5+8phi, b=59+89phi, D=89+148phi,
    Vb=5556+8989phi, Xb=62229+100689phi.

Direct multiplication gives `b Vb^2-phi Xb^2=1`. The norm of b is811;
trial division only through23 suffices since sqrt811<29 and the next
prime is29. Thus b generates a prime ideal. The inherited primitive-norm,
golden PID and positive-unit-square results leave exactly norm2 and
norm2b mirrors, after harmless unit rescaling. This use of the literal
lattice is essential; no class-wide nonreflectivity transfer is made.

At the inherited actual E8+A1(c) endpoint, the outgoing edge retains
E7+A1(c). Finite Coxeter classification leaves a norm2 opposite which
is orthogonal, joins c by3 or5, or extends E7 to E8 at node6. A connected
finite rank9 diagram cannot contain E7. For a norm2b root, reflection
integrality and the prime over811 forbid every nonzero finite angular
invariant against the retained norm2 roots. Thus that root is orthogonal
to them. This gives exactly the candidate's U0, U1, Uphi, E8 and B rows;
it classifies first SIMPLE walls, not all roots in a possible finite
endpoint reflection group.

The collision parameters, with positive denominators, are

    U0:   tau=D X/Y;
    Ueta: tau=D U/(2Y), eta=1 or phi;
    E8:   tau=2D X/Z, Z=2Y+3;
    B:    tau=phi X/V.

Their norm equations imply the displayed physical/conjugate strips
before tau_b. In particular the E8 conjugate bound has the necessary
factor1/2. The k5 parity is also correct: a=1 modulo2, so U1 has
U=1 modulo2 and Uphi has U=phi modulo2. This differs from the k4 parity.
The conservative coefficient bounds B<=8000,4000,14000,9500,10051 and
the maximum2,1,3,3,2 possible integer A values per B give exactly

    16002+4001+42003+28503+20104=110613

possible scalar pairs. These are mathematical domain bounds, not
runtime measurements. Positive physical orientation, exact collision,
parity, field-square and reflection-integrality filters are still needed.

## A new exact reduction: E8 is already in the U0 table

Every E8 candidate determines

    X0=2X, Y0=Z=2Y+3,
    Y0^2-D X0^2=1.

This is precisely the U0 equation, at precisely the SAME collision
`D X0/Y0`. Both its physical and conjugate scalar domains are the U0
domains. Conversely a U0 solution gives an integral E8 candidate exactly
when

    X0 is in2O, Y0-3 is in2O, and (Y0-3)/2>=0.

It then reconstructs `X=X0/2`, `Y=(Y0-3)/2`. Therefore a checker can
compute the U0 field squares once and reconstruct every E8 tie by these
tests. It must not delete the E8 wall as a duplicate mirror: the normals
are different even though their collision is identical.

This distinction has geometric force. For such a pair,

    m_U0=2 m_E8 + sum_{j=0}^6 t_j (-e_j),
    (t_0,...,t_6)=(2,3,4,6,5,4,3).

The U0 inequality is redundant in the E8 corner. Their product is+1,
so both cannot be distinct simple walls of an acute chamber. A routine
which chooses the orthogonal U0 row first would misidentify the full
germ despite obtaining the correct edge endpoint.

## A complete elementary exact checker

All required arithmetic is in Z[phi]; no Pell prefix or floating square
oracle is needed. For a scalar z=A+Bphi in a strip |z'|<c, put

    A0=floor(B(phi-1))=(isqrt(5 B^2)-B)//2.

For B=0 the same formula gives0. Enumerating A0-1,...,A0+2 covers all
five strips since c<=5/4; exact conjugate-sign filtering then yields the
smaller counts above. This gives constant working storage even with a
straightforward implementation.

For any required square radicand w=A+Bphi in O, use

    N(w)=A^2+AB-B^2, Tr(w)=2A+B.

If N(w) is not a nonnegative integer square, reject. Otherwise for
each possible norm s=+/-sqrt(N(w)), a square root v=a+bphi must satisfy

    (2a+b)^2=Tr(w)+2s,
    b^2=(Tr(w)-2s)/5.

Test both signs of each integral square root, require equal parity of
2a+b and b, reconstruct a, and verify v^2=w exactly. This is sufficient
as well as necessary, including the trace-zero case. Nonintegral
radicands are rejected when an O-valued root is required. The Ueta
division by4 is therefore checked rather than silently rounded. Actual
normal/reflection integrality and primitive content are verified in the
original binary-plus-native-E8 basis.

Use exact signs of quadratic elements to retain all candidates with
`0<tau<=tau_b`, select the least tau, and retain EVERY tied normal.
The finite first-simple classification proves this minimum is the actual
edge endpoint: if an unlisted root stopped the edge earlier, that finite
endpoint would have a first actual simple opposite in one of the listed
profiles, contradicting the exhaustive scalar test. This argument relies
on the inherited actual starting germ and local finiteness; it does not
infer an actual edge from a positive definite formal Gram alone.

For the full oriented endpoint germ, a finite exact redundancy check
is enough. For each possible new simple q among the ties, use the basis
of the eight retained normals plus q. Check that every tied candidate
normal is a nonnegative combination in that basis and that q has its
claimed native finite products. The actual new simple is necessarily
among the ties; these containments identify its chamber cone and expose
the U0/E8 redundancy above. A root in that cone supplies no extra facet.
Retain ambiguity or a failed containment as a genuine unresolved output;
do not silently choose a row by table order.

A certificate should bind the literal source/inputs, complete loop
domains and filters, per-type counts, all surviving scalars/square
witnesses, exact minimum/ties, original-basis normals and the final cone
containment witnesses. Discarded non-square scalars can be recomputed
by the same bounded independent checker; a huge corpus is unnecessary.
The 110613-pair bound permits a narrow stdlib worker with no eager JSON,
but no time/RSS or execution admission is asserted here.

## Conditional obstruction and the actual next decision

The candidate's map sends the orthogonal frames

    (r,g,Pold) -> (mb,c,Pnew),
    Gram=diag(2b,2,-2/phi).

Its formulas are integral on the original full lattice, not merely its
orthogonal K-form: modulo2, u=phi^2, Vb=phi, Xb=phi^2 give uVb=1 and
phi Xb=1, making both displayed halves integral. The determinant is1,
so the inverse is integral. Direct trace expansion gives

    tr(T)=7+phi b Vb/2-phi Xb-a(Xb+1)/2
         =253415+410024phi>10.

The old actual vertex follows from the inherited full Bugaenko facet
and the acute finite actual walls, not from form similarity. IF the
complete edge certificate establishes mb as first with the stated full
germ, T matches the full oriented actual germs. The integral-arrangement
principle then makes it a chamber stabilizer of infinite order. Under
the inherited arithmetic-lattice/discreteness inputs, this excludes the
literal host. Until that finite firstness gate passes it is a conditional
obstruction only. An earlier unit or B winner must be retained and its
own actual continuation studied; it is not discarded by this trace.

Thus the next substantive obligation is the single complete bounded
edge certificate above. It either closes k5 through the displayed period
or returns an actual new finite germ. Neither outcome yet constitutes a
new compact H9 example; compactness would still require a complete
all-ray/facet certificate if the route continues.
