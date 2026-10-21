# Direct full-NS maximal carriers minimalize only to \(\mathbf P^2\) or \(\mathbb F_2\)

**Status (2026-08-10): proved by the rational-surface MMP and two explicit
mirrors; no reflective-lattice classification is used.**  Let \(S\) be a
smooth projective rational surface whose complete nef cone is a hyperbolic
Coxeter chamber on the direct full Neron--Severi lattice

\[
                         L=\operatorname {NS}(S)
\]

(scaling the form by two makes no difference).  Suppose its facet group is
reflection-maximal, either among discrete reflection groups or relative to
\(O^+(L)\).  Then any minimal rational surface reached by a birational
morphism from \(S\) is

\[
                         \boxed{\mathbf P^2\quad\hbox{or}\quad\mathbb F_2.}
                                                                    \tag{0.1}
\]

Moreover, in the \(\mathbb F_2\) case no centre of the chosen blowdown
presentation lies on the successive strict transform of the negative
section.

This is a direct/full-lattice theorem.  In an adapted host containing
half-boundary classes, the ruling-swap reflection has an additional parity
test and a half-section can repair a divisibility failure.  No adapted-host
version of (0.1) is claimed.

## 1. Exact theorem

Let

\[
                         \pi:S\longrightarrow S_0          \tag{1.1}
\]

be a birational morphism to a minimal smooth rational surface.  Assume:

1. the complete Mori cone of \(S\) is rational polyhedral, so every
   irreducible negative curve supplies its literal nef facet;
2. reflection in every nef facet preserves \(L\);
3. the projectivized complete nef cone is a Coxeter fundamental chamber
   for its facet group \(W\subset O^+(L)\); and
4. \(W\) has no proper reflection overgroup in \(O^+(L)\).  Absolute
   discrete-reflection maximality, or maximality among arithmetic
   reflection groups when the chamber has finite volume, also suffices.

### Theorem 1.1 (minimal-model reduction)

Under these assumptions, \(S_0\) is \(\mathbf P^2\) or \(\mathbb F_2\).
If \(S_0=\mathbb F_2\), every centre of (1.1) avoids the successive strict
transform of the negative section.

The proof uses only the standard classification of minimal rational
surfaces as \(\mathbf P^2\) and \(\mathbb F_e\) with \(e=0\) or \(e\ge2\),
plus the following two elementary calculations.

## 2. Every blowup of \(\mathbb F_0\) has a missing mirror through Amp

Write \(c,d\) for the two ruling classes on
\(\mathbb F_0=\mathbf P^1\times\mathbf P^1\):

\[
 c^2=d^2=0,\qquad c.d=1.                                 \tag{2.1}
\]

Use the same letters for their pullbacks to \(S\), and put

\[
                              r=c-d.                       \tag{2.2}
\]

Then \(r^2=-2\), and

\[
                     s_r(x)=x+(x.r)r                      \tag{2.3}
\]

is an integral reflection of \(L\).  It fixes every total exceptional
coordinate and swaps the two ruling classes.

The mirror \(r^\perp\) meets the ample-cone interior after an arbitrary
sequence of blowups, including infinitely-near blowups.  Indeed, choose a
\(\pi\)-ample rational exceptional divisor \(D\).  For \(N\gg0\),

\[
                         A=N\pi^*(c+d)+D                   \tag{2.4}
\]

is ample.  Exceptional classes are orthogonal to \(c,d\), and
\((c+d).(c-d)=0\), so

\[
                              A.r=0.                       \tag{2.5}
\]

Because Amp\((S)\) is open, the mirror crosses its interior.  A reflection
belonging to the chamber group cannot have a mirror through the interior of
its fundamental chamber.  Hence \(s_r\notin W\), while

\[
                         \langle W,s_r\rangle\subset O^+(L)              \tag{2.6}
\]

is a discrete reflection group properly containing \(W\).  This
contradicts reflection maximality and excludes \(\mathbb F_0\), uniformly
in the number and incidence of the blowups.

This strengthens the earlier ruling-swap obstruction for one explicit
affine-tail family: neither that tail nor any special choice of centres is
used here.

## 3. The negative section excludes \(\mathbb F_e\) for \(e>2\)

Let \(s,f\) be the negative section and a fibre on \(\mathbb F_e\):

\[
                         s^2=-e,\qquad f^2=0,\qquad s.f=1. \tag{3.1}
\]

Suppose exactly \(t\) blowup centres lie successively on the current strict
transform of \(s\).  In orthogonal total-transform coordinates its final
strict transform is

\[
                    \widetilde s=s-\sum_{i\in I}e_i,
                    \qquad \widetilde s^2=-(e+t),
                    \qquad f.\widetilde s=1.              \tag{3.2}
\]

It is a primitive irreducible negative curve, hence gives a facet of the
complete nef cone.  Reflection in this facet sends

\[
                f\longmapsto
                f+{2\over e+t}\widetilde s.               \tag{3.3}
\]

Since \(\widetilde s\) is primitive, this preserves \(L\) only if

\[
                              e+t\mid2.                    \tag{3.4}
\]

Thus every \(e>2\) is impossible already from facet integrality.  When
\(e=2\), equation (3.4) forces \(t=0\), proving the final assertion of
Theorem 1.1.  The surface \(\mathbb F_1\) is not minimal: its negative
section contracts to \(\mathbf P^2\).

The classification of minimal rational surfaces now leaves precisely
\(\mathbf P^2\) and \(\mathbb F_2\), proving (0.1). \(\square\)

## 4. Why the host qualification is essential

For an adapted host

\[
 M=\left\langle L,B_1/2,\ldots,B_b/2\right\rangle,
 \qquad q_M=2q_L,                                       \tag{4.1}
\]

the vector \(r=c-d\) still defines the same reflection, but

\[
 s_r(B_\alpha/2)
 =B_\alpha/2+{B_\alpha.r\over2}r.                         \tag{4.2}
\]

Consequently

\[
 s_r\in O(M)
 \quad\Longleftrightarrow\quad
 B_\alpha.(c-d)\equiv0\pmod2
 \quad\hbox{for every }\alpha.                            \tag{4.3}
\]

Odd boundary parity can therefore remove the ruling-swap reflection from
the supplied adapted lattice.  Likewise, adjoining a half of a negative
section can change the divisibility in (3.3).  Reflection maximality never
authorizes enlarging the declared host first.  The direct theorem must not
be silently imported into the affine-H3 adapted-host lane.

There is nevertheless a sharp necessary exponent-two bound.  Every such
adapted host satisfies

\[
                            L\subseteq M\subseteq {1\over2}L.             \tag{4.4}
\]

If reflection in the primitive section \(\widetilde s\), of square
\(-n=-(e+t)\), preserves \(M\), then (3.3) belongs to \(M\).  Hence
\((2/n)\widetilde s\in(1/2)L\), and primitivity gives

\[
                              \boxed{n\mid4.}              \tag{4.5}
\]

Thus an exponent-two adapted host can permit only section squares
\(-1,-2,-4\).  This condition is necessary, not sufficient: the actual
half-class must belong to the declared boundary glue.

## 5. Consequence for the remaining vertical-depth problem

Any direct full-NS attempt at an unbounded reflection-maximal vertical-tail
family may now start only from a plane model or from \(\mathbb F_2\), with
all centres away from its negative section.  This is a genuine global MMP
reduction, but it is not yet a depth bound.  A single free proximity chain
over either surviving core can still pass the local canonical-color and
signed-permutation saturation tests.  Bounding its length requires a mixed
root or a total marked-degree argument specific to these two cores.

For affine H3, the natural host contains the classes \(B_i/2\); equation
(4.3), not (0.1), is the exact available test.  Thus the genus-three stable
invariant no-go and the adapted parity gate remain logically separate from
the direct minimal-model reduction proved here.
