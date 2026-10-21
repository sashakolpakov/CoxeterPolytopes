# Hostile audit of the first-mixed-ray standard-chain bound

## Verdict

**PASS.**  After an independent reconstruction, the conditional depth
theorem is correct and is stronger than the first direct-lattice draft.  On
the exact marked plane frame, finite volume forces \(m\leq11\); on the exact
marked \(\mathbb F_2\) frame it forces \(m\leq10\).  Since the corresponding
hyperbolic dimensions are \(m\) and \(m+1\), respectively, both lanes satisfy

\[
                              \boxed{\dim P\leq11}.
\]

The final proof is parity-free.  It applies to a direct full-Neron--Severi
host and to an adapted exponent-two host with the same exact physical frame.
No reflective-lattice classification, Nikulin--Vinberg bound, chamber census,
or surface census enters the proof.

This remains a conditional marked-frame theorem.  Reflection maximality by
itself does not force either marked frame.  This is not an existence or sharpness theorem
at dimensions ten or eleven.

## 1. The two spacelike vectors are exact

Use orthogonal total-transform coordinates on the plane blowup.  For
\(m\geq12\), the proposed vector is

\[
 w_{\mathbf P^2}=7H-3e_1-2e_2-\cdots-2e_{12}.
\]

Direct intersection gives

\[
 w_{\mathbf P^2}^2=-4,\qquad K.w_{\mathbf P^2}=4.
\]

Against the marked walls

\[
 H-e_1-e_2-e_3,\quad e_i-e_{i+1},\quad e_m,
\]

all pairings vanish except

\[
 w_{\mathbf P^2}.(e_1-e_2)=1,\qquad
 w_{\mathbf P^2}.\alpha_{12}=2.
\]

Here \(\alpha_{12}\) is the endpoint \(e_{12}\) when \(m=12\), and
\(e_{12}-e_{13}\) when \(m>12\).  Thus every displayed standard wall has
nonnegative pairing.

On \(\mathbb F_2\), write \(s^2=-2\), \(f^2=0\), and \(s.f=1\).  For
\(m\geq11\),

\[
 w_{\mathbb F_2}=4s+9f-2e_1-\cdots-2e_{11}=f+2u_{11}
\]

has

\[
 w_{\mathbb F_2}^2=-4,\qquad K.w_{\mathbb F_2}=4.
\]

Its only nonzero standard-wall pairings are

\[
 w_{\mathbb F_2}.s=1,qquad
 w_{\mathbb F_2}.\delta_{11}=2.
\]

The independent hostile tests reconstruct these vectors and all standard
pairings through rank forty without importing the producer certifier.

## 2. Every plane physical wall is nonnegative

Let a nonexceptional irreducible curve distinct from the marked line have
class

\[
 C=dH-\sum_i\mu_i e_i,qquad
 \mu_1\geq\mu_2\geq\cdots\geq0.
\]

Monotonicity follows equally from proximity or from nonintersection with the
displayed exceptional primes.  Nonintersection with the unblown line gives

\[
                    d\geq\mu_1+\mu_2+\mu_3.          \tag{2.1}
\]

Put \(q=-C^2\in\{1,2,4\}\) and
\(T=\sum_{i>12}\mu_i\).  Smooth-rational adjunction and direct pairing give

\[
 \sum_i\mu_i=3d+q-2,
 \qquad
 C.w_{\mathbf P^2}=d-\mu_1-2q+4+2T.                 \tag{2.2}
\]

For \(q=1\), the right side is at least two.  For \(q=2\), it is at least
zero.  The only serious case is a possibly primitive \((-4)\)-curve.  Set

\[
 R=d-\mu_1.
\]

Then (2.1)--(2.2) say

\[
 R\geq\mu_2+\mu_3,qquad C.w_{\mathbf P^2}=R-4+2T.  \tag{2.3}
\]

No divisibility step is needed:

- If \(T\geq2\), (2.3) is nonnegative.
- If \(T=1\), integral monotonicity forces
  \(\mu_2,\mu_3\geq1\), hence \(R\geq2\), and (2.3) is again nonnegative.
- If \(T=0\), all support lies in the first twelve places.  The square and
  adjunction equations become

  \[
  \sum_{i=2}^{12}\mu_i=2d+R+2,qquad
  \sum_{i=2}^{12}\mu_i^2=2dR-R^2+4.                 \tag{2.4}
  \]

  Negative pairing would require \(R=0,1,2\), or \(3\).  For \(R=0,1\),
  the line and monotonicity bounds make the first sum in (2.4) impossible.
  For \(R=2\), the only sum candidate is \(d=3\) with ten ones among the
  eleven remaining places; their square-sum is ten rather than the required
  twelve.  For \(R=3\), the remaining total is at most twelve, while (2.4)
  and the small values of \(d\) rule out every case.  Hence \(R\geq4\).

The line and every exceptional prime were checked separately.  This exhausts
the complete smooth-rational \((-1),(-2),(-4)\) facet alphabet and proves

\[
                         C.w_{\mathbf P^2}\geq0
\]

for every physical wall.

## 3. Every ruled physical wall is nonnegative

For a hidden irreducible strict transform write

\[
 C=as+bf-\sum_i\mu_i e_i,qquad
 \mu_1\geq\mu_2\geq\cdots\geq0.
\]

The negative section and marked strict fibre give

\[
                         b\geq2a,qquad
                         a\geq\mu_1+\mu_2.            \tag{3.1}
\]

For \(q=-C^2\in\{1,2,4\}\) and
\(T=\sum_{i>11}\mu_i\), adjunction gives

\[
 \sum_i\mu_i=2b+q-2,qquad
 C.w_{\mathbb F_2}=a-2q+4+2T.                        \tag{3.2}
\]

The \(q=1,2\) cases are nonnegative because \(a=C.f\geq0\).  For \(q=4\),
the pairing is \(a-4+2T\).  If \(T\geq2\), it is nonnegative.  If \(T=1\),
monotonicity and (3.1) give \(a\geq2\), again enough.  When \(T=0\),

\[
 \sum_{i=1}^{11}\mu_i=2b+2,qquad
 \sum_{i=1}^{11}\mu_i^2=-2a^2+2ab+4.               \tag{3.3}
\]

The four negative-pairing possibilities \(a=0,1,2,3\) all fail.  The first
two contradict the sum and \(b\geq2a\).  For \(a=2\), only \(b=4\) and ten
ones survive the sum bound, but their square-sum is ten instead of twelve.
For \(a=3\), the largest allowed monotone pattern has total twelve, whereas
\(b\geq6\) requires at least fourteen.  Thus \(a\geq4\).

The section, marked fibre, and exceptional primes have already been checked,
so every ruled physical facet is nonnegative on \(w_{\mathbb F_2}\).

## 4. The wall-cone and finite-volume implication is sound

Let \(\mathcal D\) be the homogeneous cone cut out by all oriented physical
facet inequalities.  An ample class \(h\) lies in
\(\operatorname{int}\mathcal D\).  Either mixed vector \(w\) above lies in
\(\mathcal D\) and has negative square.  Since

\[
 h^2>0,\qquad w^2<0,
\]

the path \(h+tw\) has a first positive parameter \(t_0\) for which

\[
                         z=h+t_0w,qquad z^2=0.
\]

For any convex cone,

\[
 \operatorname{int}\mathcal D+\mathcal D
 \subseteq\operatorname{int}\mathcal D.
\]

Consequently \(z\) is an interior null direction of the wall cone, not an
ordinary isolated cusp.  A projective neighborhood of \(z\) contains an
open ideal-boundary patch.  The hyperbolic radial integral over such a patch
diverges, so the chamber has infinite volume.

There is also a useful AG cross-check.  If the declared homogeneous cone is
literally the complete nef cone, then \(w\in\operatorname{Nef}(S)\) would
contradict Hodge immediately, because every nef class has nonnegative square.
The wall-cone proof is the weaker formulation and does not incorrectly call
\(w\) itself nef before cone completeness has been invoked.

Scaling the form or passing to an adapted real host changes neither the wall
hyperplanes nor the sign of \(w^2\).  This proves that the final theorem is
genuinely parity-free and covers adapted hosts on the same exact frames.

## 5. Scope failures were tested rather than hidden

The hypotheses used above are real geometric content.

1.  The plane class \(\alpha_0=H-e_1-e_2-e_3\) must be an actual unblown
    irreducible line wall.
2.  The exceptional primes must form the displayed free chain, so their
    nonintersection inequalities impose monotone multiplicities.
3.  On \(\mathbb F_2\), the second point is the fibre-following free infinitely-near
    point where the strict fibre meets \(E_1\).  Only \(E_1\)
    is exceptional there, so it is not satellite.  No centre lies on the
    negative section.
4.  Every hyperbolic facet must belong to the complete physical
    \((-1),(-2),(-4)\) alphabet.  Other negative squares, higher-genus walls,
    or a weak incomplete divisor subsystem are outside the theorem.

Independent formal countermodels show the necessity of these inequalities.
For example,

\[
                 2H-e_1-\cdots-e_8
\]

has square \(-4\), canonical pairing \(2\), and negative pairing with
\(w_{\mathbf P^2}\), but it intersects the marked line negatively.  Moving
one of its simple multiplicities past a zero repairs the line pairing but
violates free-chain monotonicity.  Similarly,

\[
                 s+2f-e_1-\cdots-e_6
\]

has the correct rational \((-4)\) numerics and negative ruled mixed pairing,
but intersects the marked fibre wall negatively.  A reordered six-point
version repairs both global-wall pairings only by violating proximity.
These are numerical scope regressions, not claims of simultaneous effective
geometric realizations.

Thus the theorem gives a genuine AG dimension bound on two exact lanes, but
it does not prove a universal rational-carrier bound.  The remaining bridge
is to force one of the marked line/fibre configurations, or an equally strong
replacement inequality, from an arbitrary one-chain maximal carrier.

## 6. Independent replay and frozen producer bytes

The hostile regression is
[test_rational_single_chain_first_mixed_ray_bound_hostile.py](../tests/test_rational_single_chain_first_mixed_ray_bound_hostile.py).
It imports no producer helper.  It independently reconstructs both Lorentz
forms, both mixed rays, every standard wall pairing, the complete finite
low-\(R\) and low-\(a\) cases, and the formal scope countermodels.

The audited producer bytes are frozen at:

```text
results/rational_single_chain_first_mixed_ray_bound.md
  57bef687ef27540c92fe277e722200f6906391a2584dbd62c8a680e1367bdbd1
scripts/certify_rational_single_chain_first_mixed_ray_bound.py
  f9c4c4bcaa9dfb0d758daff31e6ae0ad6399614d11b047cd8ca571c834564abb
data/enumerations/rational_single_chain_first_mixed_ray_bound.json
  25e0ba8392719b68e4e7e0d328bd6795419c48de620d7e84e796749e2acee078
tests/test_rational_single_chain_first_mixed_ray_bound.py
  db8e9edcb8c8af5edb44c0ed11f98ec37d93c86ddf58e0f22191d6d5559c435c
```
