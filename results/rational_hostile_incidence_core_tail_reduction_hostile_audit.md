# Hostile audit of the rational hostile-incidence synthesis

**Verdict (2026-08-10): PASS after four scope repairs.** The synthesis now
faithfully composes the frozen weak-incidence, finite-vertex, minimal-H1,
minimal-H2, and H3 reductions. It proves no universal rank theorem. In
particular, it neither turns a bounded relative core into a bound on the
original Picard rank nor turns the elementary finite-volume span condition
into a bound on a signed-permutation tail.

The repairs were:

1. the finite-vertex statement now retains the SNC, negative-definite,
   normal projective, and Q-Gorenstein contraction hypotheses;
2. H2 discrepancy language now distinguishes Grauert's analytic contraction
   from a supplied Q-Gorenstein contraction;
3. the spherical H3 discrepancy vector is conditional on a supplied
   Q-Gorenstein contraction, while negative definiteness alone supplies only
   the complex-analytic contraction;
4. affine H3 semiampleness is explicitly separate from the global tail
   lemma, and the finite-volume paragraph now warns that one mixed
   W(B_m)-orbit already spans an arbitrary tail.

No component theorem was strengthened during synthesis.

## 1. Dependency boundary

The synthesis uses five independently audited components.

| component | exact imported conclusion | conclusion not imported |
|---|---|---|
| weak wrong-log lane | under complete Mori generation and W1--W3, ρ≤26 and n≤25 | raw essential surjectivity |
| finite-vertex lane | an SNC negative packet with the declared normal projective Q-Gorenstein contraction is ε-klt for every ε<1/2 | rank or tail length |
| minimal H1 | for B₁²=B₂²=-4 and B₁·B₂=2 with two transverse contacts, ρ_core≤14, sharply | higher/parallel/dotted incidence or original rank |
| minimal H2 | for B²=-4, R²=-2, and B·R=2, ρ_core≤12, sharply | higher/dotted incidence, a complete equality carrier, or original rank |
| H3 outside H1/H2 | the exact Q=3, Q=4, Q>4 star trichotomy | affine semiampleness or a tail bound |

The synthesis table and terminal section now preserve all five boundaries.

## 2. H1 and H2 are only the minimal finite-angle cases

In the standard doubled host, a square-minus-four physical curve has root
h=B/2 of norm -2.

For H1, two such roots have host pairing

    q_M(B₁/2,B₂/2) = (B₁·B₂)/2.

Standard-host integrality makes the physical intersection I even. The case
I=2 has normalized product 1/2 and Coxeter order three. The case I=4 is
parallel, and I>4 is ultraparallel/dotted. The genus-one adjoint proof uses
exactly I=2: for I>2 the dualizing degree I-2 is positive on each component.
The synthesis correctly claims no higher-incidence H1 theorem.

For H2, the roots B/2 and R have norms -2 and -4 and host pairing B·R=I.
The value I=2 has normalized product 1/sqrt(2) and Coxeter order four.
The value I=1 is not a Coxeter angle in this standard crystallographic host,
and I≥3 is ultraparallel. Again, only I=2 enters the rank-12 adjoint core.
The cohomological H2 core permits a nonreduced length-two intersection, but
the lc singularity statement requires two transverse contacts; an ordinary
tangency is non-lc after resolving it. The synthesis keeps that distinction.

The numerical core identities replay exactly:

    H1: (K+B₁+B₂)² = 14-ρ_core,
    H2: (K+B+R)²   = 12-ρ_core.

The free exceptional complement remains an orthogonal <-2>^m host summand,
so these are core bounds only.

## 3. H3 contraction and abundance qualifications

Put t_i=B_i·E, T=Σt_i, Q=Σt_i², and

    F = 4E + Σt_i B_i.

Under the H3 hypotheses, direct physical intersection gives

    F·B_i=0,   F·E=Q-4,   F²=4(Q-4),   K·F=-4+2T.

Completeness of the Mori alphabet is essential to infer that F is nef when
Q≥4.

For Q=3, the profile is (1,1,1). The physical packet is negative definite
and its fundamental cycle 3E+B₁+B₂+B₃ has arithmetic genus one. Grauert
gives a complex-analytic contraction over C. The fundamental-cycle
calculation makes any supplied normal contraction non-rational; only if it
is Q-Gorenstein may one solve for discrepancies

    (-2,-1,-1,-1),

and conclude that it is not log canonical. The repaired synthesis says
precisely this.

For Q=4, H3 forces the four-unit profile. The nef isotropic class F is not
known semiample. Conditional on F=mG defining a fibration, adjunction gives
only (m,g)=(1,3) or (2,2), and relative minimalization gives core ranks at
most 26 or 18. The vertical blowup count remains unbounded. Abundance and
tail length are two separate gates.

For Q>4, the H3 component theorem independently supplies a projective
contraction of the disjoint B_i using a finite-Mori-generator
base-point-free argument; it does not reuse the failed weak-lane nef class.
On the quotient,

    Ebar²=Q/4-1>0,   K_Y·Ebar=T/2-1>0.

The supported-boundary inequality and the wrong Hodge-index direction are
therefore faithfully imported. No orbifold-c2 rank bound is asserted.

## 4. Canonical coloring and every factor of two

The standard host convention is q_M=2q_NS. Applying it to the physical
wall root, not to the physical curve indiscriminately, gives:

| physical curve | host root | q_M(root,root) | q_M(K,root) |
|---|---|---:|---:|
| R²=-2 | R | -4 | 0 |
| E²=-1 | E | -2 | -2 |
| B²=-4 | B/2 | -2 | +2 |

Thus only the norm-minus-four physical roots are canonically color zero,
and they span a subspace of K-perp of rank at most ρ-1. The synthesis does
not confuse B with B/2.

The H3 host vector is also scaled correctly:

    e=E,   h_i=B_i/2,
    r=F/2=2e+Σt_i h_i,
    r²=2(Q-4),   r·e=Q-4,   r·h_i=0.

For a free m-step exceptional tail, the norm-minus-four roots
e_i-e_{i+1} span only the sum-zero A_{m-1} hyperplane. In total-transform
coordinates the canonical tail component is e₁+...+e_m, and the terminal
norm-minus-two root e_m supplies the missing facet direction. This is why
canonical coloring is a genuine obstruction for some full chambers but
does not eliminate the B_m tail.

## 5. What finite volume forces—and what it does not

Suppose the carrier lattice split orthogonally as a Lorentzian core plus a
negative m-dimensional tail, and every reflection root lay wholly in one
factor. The generated group would preserve the corresponding proper
totally geodesic core. In Fermi coordinates around it, hyperbolic volume has
a radial factor of the form

    cosh(r)^a sinh(r)^(m-1) dr.

Every nonzero spherical tail sector therefore has infinite radial volume.
A finite-covolume chamber cannot be a pure orthogonal product: some root or
group element must mix core and tail.

That observation is much weaker than a tail bound. The local tail group is
the full signed-permutation W(B_m). If v is any nonzero tail projection,
choose a nonzero coordinate v_i. The difference between v and its
single-coordinate sign change is 2v_i e_i; permutations then put every
coordinate vector in the orbit span. Hence one mixed-root orbit already
spans the entire tail for every m. Finite volume supplies an orbit-span
condition, not a bound on m, on mixed-root coefficients, or on the
adjunction/effectivity complexity of those roots.

The synthesis now makes exactly this claim. Its “global mixed-root tail
lemma” remains a requested theorem, not a proved consequence of finite
volume or reflection maximality.

## 6. Final scope

The audited synthesis establishes:

- the complete weak lane has n≤25;
- the minimal H1 and H2 hostile cores have sharp ranks 14 and 12;
- H3 has its exact spherical/affine/hyperbolic split;
- pure exceptional reflection saturation cannot bound the erased tail; and
- after the separate affine abundance gate, the remaining common problem is
  a global geometric bound on mixed core/tail attachment.

It does not establish a universal dimension bound, higher-incidence H1/H2,
affine H3 abundance, trace descent, a uniform total facet-degree bound, or
an unbounded reflection-maximal counterfamily.

The dependency hashes and the synthesis scope are pinned independently by
tests/test_rational_hostile_incidence_core_tail_reduction_hostile.py. Those
tests recompute the finite-angle normalizations, both core constants, the H3
star identities, the canonical-color factors, and the B_m orbit-span
warning without importing any component certifier.
