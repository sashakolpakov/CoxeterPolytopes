# The first mixed ray bounds the two standard minimal-base chains

**Status (2026-08-10): exact conditional depth theorem; hostile review
pending.**  On the two marked minimal-base frames isolated by the rational
surface reduction, the first mixed hyperideal ray cannot be cut by any
physical wall.  Therefore

\[
 \boxed{m\le11\quad\text{over }\mathbf P^2,\qquad
        m\le10\quad\text{over }\mathbb F_2.}           \tag{0.1}
\]

Both bounds give hyperbolic dimension at most eleven.  The proof uses only
intersection theory, proximity, smooth-rational adjunction, and the
elementary finite-volume cone criterion.  It uses no
Nikulin--Vinberg classification and no chamber or surface census.

The theorem is deliberately about the **exact marked standard frames**.  It
does not prove that reflection maximality forces an arbitrary one-chain
blowup into either frame.  Unlike the preceding doubled-wall argument, the
proof below uses no parity or divisibility: it applies on either the direct
full-NS host or an adapted half-boundary host once the same physical frame
and complete alphabet have been constructed.

## 1. Hypotheses and the finite-volume cone lemma

Let \(S\) be a smooth projective rational surface whose projectivized
complete nef cone is a finite-sided finite-volume Coxeter chamber.  Assume
that every hyperbolic nef facet is represented by a smooth rational curve
of square \(-1,-2\), or \(-4\).  The declared reflection host may be the
direct Neron--Severi lattice or an adapted exponent-two host.  The proof
uses the curves themselves, not their primitive wall roots.

We use the following standard projective criterion.

### Lemma 1.1 (spacelike wall-cone obstruction)

Let \(P\) be a finite-sided hyperbolic polyhedron and orient its facet
normals \(r_i\) so that the homogeneous inequality cone

\[
 \mathcal C(P)=\{x:x.r_i\ge0\text{ for every }i\}      \tag{1.1}
\]

meets the chosen positive sheet in the cone over \(P\).  If a spacelike
vector lies in \(\mathcal C(P)\), then \(P\) has infinite volume.

### Proof

Choose a future timelike \(h\) in the interior of \(\mathcal C(P)\).  If a
spacelike \(w\in\mathcal C(P)\) existed, the ray \(h+tw\) would first meet
the future null cone at some \(z=h+t_0w\), \(t_0>0\).  For every convex
cone,

\[
             \operatorname {int}\mathcal C(P)+\mathcal C(P)
             \subset\operatorname {int}\mathcal C(P),
\]

so \(z\) is an **interior** null direction of the wall cone.  An open patch
of the ideal boundary consequently lies inside the chamber's projective
closure.  In the ball model the hyperbolic density has a nonintegrable pole
on such a positive-measure boundary patch, so the volume is infinite.  This
contradiction proves the lemma. \(\square\)

This is the same homogeneous-cone criterion used in the preceding
standard-simplex calculation.  Notice that it concerns the complete
facet cone, not whether a spacelike test vector is itself a nef divisor.

## 2. The exact plane frame

Let \(S_m\to\mathbf P^2\) be the following marked cluster.  The first three
centres follow a line, the fourth centre leaves its strict transform, and
every exceptional successor is free on the newest exceptional component.
Assume the actual irreducible facet curves include

\[
 \alpha_0=H-e_1-e_2-e_3,
 \qquad
 \alpha_i=e_i-e_{i+1}\ (1\le i<m),
 \qquad
 \alpha_m=e_m.                                        \tag{2.1}
\]

Thus \(\alpha_0\) is an unblown strict line, all displayed exceptional
differences are the strict exceptional primes, and every other irreducible
negative curve is the strict transform of a plane curve

\[
 C=dH-\sum_{i=1}^m\mu_i e_i,
 \qquad
 \mu_1\ge\mu_2\ge\cdots\ge\mu_m\ge0.                 \tag{2.2}
\]

For \(C\ne\alpha_0\), nonintersection of distinct irreducible curves gives

\[
                    C.\alpha_0=d-\mu_1-\mu_2-\mu_3\ge0.                \tag{2.3}
\]

These proximity and line inequalities are precisely why the marked-frame
hypothesis cannot be dropped.

### 2.1 The first mixed ray

For every \(m\ge12\), put

\[
 w_{\mathbf P^2}
   =7H-3e_1-2e_2-\cdots-2e_{12}.                       \tag{2.4}
\]

This is the mixed ray \(v_1+2v_{12}\) from the standard-simplex gate.  A
direct calculation gives

\[
 w_{\mathbf P^2}^2=-4,
 \qquad K_S.w_{\mathbf P^2}=4.                        \tag{2.5}
\]

Its only nonzero pairings with (2.1) are

\[
                 w_{\mathbf P^2}.\alpha_1=1,
 \qquad          w_{\mathbf P^2}.\alpha_{12}=2.       \tag{2.6}
\]

All standard wall inequalities are therefore nonnegative.

### 2.2 Every physical hidden wall is also nonnegative

Let \(C\) in (2.2) be smooth rational and put

\[
                             q=-C^2\in\{1,2,4\},
 \qquad T_{12}=\sum_{i>12}\mu_i.                       \tag{2.7}
\]

Adjunction says

\[
 -3d+\sum_i\mu_i=q-2,
 \qquad
 \sum_i\mu_i=3d+q-2.                                  \tag{2.8}
\]

Substitution into (2.4) gives the exact identity

\[
 \boxed{
 C.w_{\mathbf P^2}
   =d-\mu_1-2q+4+2T_{12}.}                             \tag{2.9}
\]

By (2.3), \(d\ge\mu_1\).  Hence (2.9) is at least \(2\)
for a \((-1)\)-curve and at least zero for a \((-2)\)-curve.

It remains to treat \(C^2=-4\).  Put

\[
 R=d-\mu _1,\qquad T=T_{12}.
\]

Then (2.3) and (2.9) give

\[
 R\ge\mu _2+\mu _3,\qquad
 C.w_{\mathbf P^2}=R-4+2T.                            \tag{2.10}
\]

If \(T\ge2\), the pairing is nonnegative.  If \(T=1\), monotonicity
forces \(\mu _2,\mu _3\ge1\), so \(R\ge2\) and the pairing is again
nonnegative.  Suppose finally that \(T=0\).  The curve is then supported
on the first twelve exceptional places, and its square and adjunction give

\[
 \sum_{i=2}^{12}\mu_i=2d+R+2,
 \qquad
 \sum_{i=2}^{12}\mu_i^2=2dR-R^2+4.                   \tag{2.11}
\]

The four possible negative-pairing values \(R=0,1,2,3\) are impossible:

- If \(R=0\), the line inequality forces every \(\mu_i\) with \(i\ge2\)
  to vanish, contradicting the first equation in (2.11).
- If \(R=1\), monotonicity and \(\mu_2+\mu_3\le1\) give total remaining
  multiplicity at most one, again contradicting (2.11).
- If \(R=2\), either \(\mu_3=0\), which leaves total remaining
  multiplicity at most two, or \(\mu_2=\mu_3=1\) and every remaining
  multiplicity is at most one.  The latter case forces \(d=3\) and ten
  ones among the eleven remaining places.  Their square-sum is ten, while
  the second equation in (2.11) requires twelve.
- If \(R=3\), the eleven remaining places have total multiplicity at most
  twelve.  For \(d=3\), \(\mu_1=0\) forces them all to vanish; for \(d=4\),
  \(\mu_1=1\) lowers their maximum total to eleven; and for \(d\ge5\),
  (2.11) requires \(2d+5\ge15\), already above twelve.

Thus \(R\ge4\), so every physical \((-4)\)-wall also has nonnegative
pairing with \(w_{\mathbf P^2}\).  Notice that this argument never divides
the curve class by two.

The exceptional primes and the line are already the walls (2.1), so this
exhausts the complete physical alphabet.  Equations (2.5)--(2.11) put a
spacelike vector in the complete homogeneous wall cone.  Lemma 1.1 proves
that finite volume is impossible for \(m\ge12\).  Therefore

\[
                              \boxed{m\le11.}           \tag{2.15}
\]

## 3. The exact \(\mathbb F_2\) frame

Let \(s^2=-2\), \(f^2=0\), and \(s.f=1\).  Blow up a point on a fibre, then
blow up the fibre-following free infinitely-near point where its strict
transform meets the first exceptional curve.  Only \(E_1\) is exceptional
at that second centre, so the point is not satellite.  Every later centre
is free on the newest exceptional component and every centre is off \(s\).
Assume the actual irreducible
facet curves include

\[
 s,qquad
 \delta_0=f-e_1-e_2,qquad
 \delta_i=e_i-e_{i+1}\ (1\le i<m),qquad
 \delta_m=e_m.                                        \tag{3.1}
\]

The fibre-following qualification makes both \(\delta_0\) and
\(\delta_1=e_1-e_2\) physical.  The subsequent free-tail qualification
makes every later exceptional difference physical.

Every hidden irreducible negative curve is a strict transform

\[
 C=as+bf-\sum_{i=1}^m\mu_i e_i,
 \qquad \mu_1\ge\mu_2\ge\cdots\ge0.                   \tag{3.2}
\]

Distinct-curve nonintersection with \(\delta_0\) and \(s\) gives

\[
              a\ge\mu_1+\mu_2,qquad b\ge2a.          \tag{3.3}
\]

The displayed section, fibre, and exceptional primes are treated directly
by (3.1).

### 3.1 The first ruled mixed ray

For every \(m\ge11\), put

\[
 w_{\mathbb F_2}
   =4s+9f-2e_1-\cdots-2e_{11}.                         \tag{3.4}
\]

This is \(f+2u_{11}\).  One has

\[
 w_{\mathbb F_2}^2=-4,
 \qquad K_S.w_{\mathbb F_2}=4.                         \tag{3.5}
\]

Its only nonzero standard-wall pairings are

\[
                       w_{\mathbb F_2}.s=1,
 \qquad               w_{\mathbb F_2}.\delta_{11}=2. \tag{3.6}
\]

### 3.2 Every physical hidden wall is nonnegative

Again put \(q=-C^2\in\{1,2,4\}\), and let

\[
                             T_{11}=\sum_{i>11}\mu_i. \tag{3.7}
\]

Since

\[
 K_S=-2s-4f+\sum_i e_i,
\]

adjunction gives

\[
                    \sum_i\mu_i=2b+q-2.               \tag{3.8}
\]

Substitution into (3.4) yields

\[
 \boxed{
 C.w_{\mathbb F_2}=a-2q+4+2T_{11}.}                   \tag{3.9}
\]

This is nonnegative for \(q=1,2\), since \(a=C.f\ge0\).

For \(q=4\), equation (3.9) reads

\[
 C.w_{\mathbb F_2}=a-4+2T_{11}.                       \tag{3.10}
\]

If \(T_{11}\ge2\), this is nonnegative.  If \(T_{11}=1\), monotonicity
forces \(\mu_1,\mu_2\ge1\); the fibre inequality in (3.3) then gives
\(a\ge2\), so (3.10) is again nonnegative.  Suppose \(T_{11}=0\).  The
curve is supported on the first eleven exceptional places, and its square
and adjunction give

\[
 \sum_{i=1}^{11}\mu_i=2b+2,
 \qquad
 \sum_{i=1}^{11}\mu_i^2=-2a^2+2ab+4.                 \tag{3.11}
\]

The four possible negative-pairing values \(a=0,1,2,3\) are impossible:

- If \(a=0\), the fibre inequality forces all multiplicities to vanish,
  contradicting \(2b+2>0\) and \(b\ge0\).
- If \(a=1\), the total multiplicity is at most one, whereas \(b\ge2\)
  makes it at least six.
- If \(a=2\), the total multiplicity is at most eleven.  Hence \(b=4\),
  and total ten forces ten ones.  Their square-sum is ten, but (3.11)
  requires twelve.
- If \(a=3\), the largest monotone pattern allowed by
  \(\mu_1+\mu_2\le3\) has total twelve, while \(b\ge6\) requires total at
  least fourteen.

Therefore \(a\ge4\), and every physical \((-4)\)-wall has nonnegative
pairing with \(w_{\mathbb F_2}\), without a divisibility assumption.

Together with (3.6), this checks every physical facet.  The spacelike vector
\(w_{\mathbb F_2}\) lies in the complete homogeneous wall cone, so Lemma
1.1 excludes finite volume for \(m\ge11\).  Hence

\[
                              \boxed{m\le10.}           \tag{3.14}
\]

## 4. Consequences and exact boundary

The Picard ranks and hyperbolic dimensions in the two frames are

\[
\begin{array}{c|c|c|c}
\text{base}&\rho(S)&\dim\mathbf H(\operatorname {NS}(S))&\text{bound}\\ \hline
\mathbf P^2&m+1&m&m\le11\\
\mathbb F_2&m+2&m+1&m\le10.
\end{array}                                             \tag{4.1}
\]

Thus both exact standard lanes satisfy the same dimension ceiling

\[
                              \boxed{\dim P\le11.}      \tag{4.2}
\]

This is a real AG depth theorem, not a census: it quantifies over every
hidden physical wall at once.  Reflection maximality is not needed after
the exact physical-frame and complete-alphabet hypotheses are imposed;
the theorem therefore applies a fortiori to a reflection-maximal chamber in
either exact frame.

The remaining global gap is now narrower and precise.

- An arbitrary one-chain blowdown need not contain the marked line wall
  \(\alpha_0\) or the marked fibre wall \(\delta_0\).
- Other proximity patterns change the monotonicity inequalities and the
  exceptional primes.
- Curves outside the complete smooth-rational \((-1),(-2),(-4)\) facet
  alphabet are not covered.

Accordingly (0.1) is not yet the universal hostile-tail theorem.  It does,
however, prove that the first mixed ray is the missing total invariant on
both exact standard frames.  Any surviving counterexample must erase one
of those global line/fibre walls or change the proximity type.  Passing to
an adapted half-boundary host alone does not evade the theorem on these
exact frames.

## 5. Exact replay

- [symbolic certifier](../scripts/certify_rational_single_chain_first_mixed_ray_bound.py)
- [JSON archive](../data/enumerations/rational_single_chain_first_mixed_ray_bound.json)
- [regression tests](../tests/test_rational_single_chain_first_mixed_ray_bound.py)

The certifier reconstructs both mixed rays through rank 32, checks their
norms and every standard-wall pairing, and freezes the parity-free
adjunction/proximity case splits species by species.  It imports no existing root
enumerator.
