# Hostile audit of the direct full-NS minimal-model reduction

## Verdict

**PASS.**  The theorem in the
[minimal-model note](rational_direct_full_ns_minimal_model_reduction.md) is
correct under its stated complete-nef and reflection-maximality hypotheses.
The rational MMP leaves only \(\mathbf P^2\), \(\mathbb F_0\), and
\(\mathbb F_e\) for \(e\geq2\).  A ruling-swap mirror excludes every
blowup of \(\mathbb F_0\), while literal facet integrality excludes
\(\mathbb F_e\) for \(e>2\) and excludes every centre on the strict
negative section of \(\mathbb F_2\).

This is a direct full-Neron--Severi result.  The adapted-host statements are
only an exact parity test for the ruling swap and the necessary condition
\(e+t\mid4\) for the section reflection.  Neither is promoted to an
adapted-host minimal-model theorem or a chain-depth bound.

## 1. The rational MMP reduction has no missing minimal base

For a smooth projective rational surface over the declared complex ground
field, a smooth minimal model is \(\mathbf P^2\) or a Hirzebruch surface
\(\mathbb F_e\), with \(e=0\) or \(e\geq2\).  The case \(e=1\) is not
minimal because its section contracts to \(\mathbf P^2\).  Thus the two
lattice obstructions in the producer note exhaust the nonplane minimal
bases; no reflective-lattice classification enters this step.

The use of a negative curve as a wall is also explicit rather than hidden:
the theorem assumes the complete rational-polyhedral Mori cone and literal
reflection integrality for every nef facet.  The strict negative section is
irreducible and negative throughout the relevant blowup presentation, so
its ray supplies exactly the facet used in the proof.

## 2. The \(\mathbb F_0\) mirror survives arbitrary blowups

In total-transform coordinates,

\[
 c^2=d^2=0,\qquad c.d=1,\qquad r=c-d,\qquad r^2=-2.
\]

Reflection in \(r\) swaps \(c,d\) and fixes every total exceptional
coordinate.  It is therefore an integral involutive isometry of the full
Neron--Severi lattice in every exceptional rank.

The geometric ample-mirror argument does not assume general-position
centres.  For an arbitrary composition \(\pi:S\to\mathbb F_0\), choose a
rational divisor \(D\), supported on the exceptional locus, which is
\(\pi\)-ample.  Relative ampleness implies that

\[
                    N\pi^*(c+d)+D
\]

is ample for sufficiently large \(N\).  Both summands pair trivially with
\(r\), so \(r^\perp\) really meets the ample interior even for satellite
and infinitely-near centres.

Since a chamber-group mirror cannot cross the interior of its fundamental
chamber, \(s_r\notin W\).  But \(s_r\in O^+(L)\), and \(O^+(L)\) is
discrete, so \(\langle W,s_r\rangle\) is a proper discrete reflection
overgroup.  This also defeats the arithmetic finite-covolume formulation.

## 3. The strict-section divisibility is exact

If \(t\) centres lie on the successive strict transform of the negative
section of \(\mathbb F_e\), then in orthogonal total-transform coordinates

\[
 \widetilde s=s-\sum_{i\in I}e_i,\qquad
 \widetilde s^2=-(e+t),\qquad f.\widetilde s=1.
\]

The coefficient of \(s\) is one, so \(\widetilde s\) is primitive even
when some of those centres are satellite points.  Put \(n=e+t\).  Its
reflection sends

\[
                       f\longmapsto f+{2\over n}\widetilde s.
\]

Primitivity shows that this vector belongs to the full Neron--Severi
lattice exactly when \(n\mid2\).  Consequently \(e>2\) is impossible, and
for \(e=2\) one must have \(t=0\).  The independent tests reconstruct the
whole rational reflection matrix and verify the if-and-only-if condition,
rather than checking only sampled output from the producer certifier.

## 4. The adapted-host qualifications are sharp

For the doubled exponent-two host

\[
 M=\langle L,B_1/2,\ldots,B_b/2\rangle,
\]

the ruling-swap shift is

\[
 s_r(B_\alpha/2)-B_\alpha/2={B_\alpha.r\over2}r.
\]

Thus all boundary pairings must be even.  Conversely, even pairings make
the reflection preserve every displayed generator.  If a pairing were odd
and \(r/2\) nevertheless lay in the integral doubled host, then
\(q_M(B_\alpha/2,r/2)=B_\alpha.r/2\) would not be integral.  Hence parity
is exact for this reflection on a legitimate integral adapted host.

For the strict section, preservation of any host satisfying
\(L\subset M\subset\tfrac12L\) implies

\[
                 {2\over n}\widetilde s\in M\subset\tfrac12L.
\]

Primitivity gives \(n\mid4\).  This is only necessary: it does not say the
required half-class occurs in the declared glue.  Odd ruling parity and the
section-glue condition therefore block import of the direct theorem into
the affine-H3 lane.

## 5. Frozen producer artifacts

The hostile replay pins:

    de94807a390d6d8ec7238b16b00deb304077cd942c757aa96dc1cba529b32abc  results/rational_direct_full_ns_minimal_model_reduction.md
    d70f4d7df8a08aa7ee219b3e71eb0382269896aad55eaf055d6faacf770f444e  scripts/certify_rational_direct_full_ns_minimal_model_reduction.py
    305c8a9a9b0bdef6ac6de34e649a7048a97e2a3b3b142c4f2210cb34a4766d87  tests/test_rational_direct_full_ns_minimal_model_reduction.py

No depth bound, adapted-host minimal-base classification, or external
reflection classification is asserted by this audit.
