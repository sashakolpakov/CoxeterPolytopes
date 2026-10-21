# Hostile audit of the affine \(B_m\) global-completion obstruction

## Verdict

**PASS after scope and certificate repairs.**  The uniform
statement in the
[producer theorem](rational_affine_b_tail_global_completion_obstruction.md)
is correct:

> For every \(m\ge2\), the explicit iterated blowup of
> \(\mathbf P^1\times\mathbf P^1\) cannot have its complete nef cone as an
> absolutely reflection-maximal Coxeter chamber on the full
> Neron--Severi host.

The proof is independent of the choice of the free centres and of any
unknown negative curves.  The class \(c-d\) is an integral reflecting root
whose mirror contains an ample class, so its reflection is absent from the
nef facet group but can always be adjoined as a proper discrete reflection
overgroup.  No Nikulin--Vinberg classification, chamber census, or
reflective-lattice list enters.

The visible fibre--section simplex also has finite volume exactly for
\(m\le4\), and the first hidden-wall inequality is correct.

The audit repaired four points which did not change those conclusions:

1. the missing TeX command in \(K.T_5=-1,\quad p_a(T_5)=0\);
2. \(T_5\) remains an integral square-minus-one class after every later
   blowup, so \(R_1\) is only the first **square-minus-two,
   minimal-height** candidate at \(m=6\), not the first integral candidate
   without qualification;
3. lattice mirrors and effective divisor classes were separated explicitly
   from irreducible smooth curves and actual nef facets.  Reflection
   maximality does not turn every lattice root into a physical wall;
4. the codimension-\(a\) jet claim was restricted to irreducible
   exact-multiplicity representatives.  The row-dependency count applies to
   every effective divisor, but exceptional fixed components invalidate the
   naive finite-branch argument for arbitrary reducible divisors.

The producer certificate and archive were regenerated so that \(T_5\) is
retained for every sampled \(m\ge5\).

The independent audit artifacts are:

- [hostile script](../scripts/audit_rational_affine_b_tail_global_completion_obstruction_hostile.py);
- [hostile archive](../data/enumerations/rational_affine_b_tail_global_completion_obstruction_hostile_audit.json);
- [hostile tests](../tests/test_rational_affine_b_tail_global_completion_obstruction_hostile.py).

The audit script imports none of the producer's functions.

## 1. Blowup geometry and the visible integral basis

In total-transform coordinates,

\[
 L=U\perp\langle-1\rangle^m,\qquad
 c^2=d^2=0,\quad c.d=1,\quad e_i^2=-1.                 \tag{1.1}
\]

The first blowup lies on the fibre \(c\).  The second centre is the point
of the first exceptional curve representing the tangent direction of
\(c\).  The other ruling fibre through the first centre has a different
tangent direction, so its strict transform

\[
                              D=d-e_1                  \tag{1.2}
\]

misses the second centre.  Every later centre lies on the newest
exceptional curve away from all intersection points, so it also misses
\(D\).  Hence \(D\) is an actual smooth rational \((-1)\)-curve for every
allowed cluster.

The other universal curves are

\[
\delta_0=c-e_1-e_2,\qquad
\delta_i=e_i-e_{i+1}\ (1\le i<m),\qquad
\delta_m=e_m.                                         \tag{1.3}
\]

They are respectively the strict transform of \(c\), the successive
exceptional strict transforms, and the last exceptional curve.  Thus they
are actual smooth rational curves, not merely numerical classes.  Direct
substitution gives

\[
\delta_0+\delta_1+2\delta_2+\cdots+2\delta_m=c          \tag{1.4}
\]

and

\[
D.\delta_1=1,                                          \tag{1.5}
\]

in addition to the affine \(B_m\) incidences.  Conversely,

\[
e_i=\delta_i+\cdots+\delta_m,\qquad d=D+e_1,
\]

and (1.4) recovers \(c\).  Therefore
\((D,\delta_0,\ldots,\delta_m)\) is an integral basis of \(L\).  The hostile
certificate checks its coordinate determinant is \(\pm1\) for every
\(2\le m\le48\), reconstructs every intersection, and verifies the Coxeter
orders \(3\) and \(4\) directly from normalized pairings.

## 2. The visible simplex and the exact \(m=4\) threshold

For the wall order

\[
                        (D,\delta_0,\delta_1,\ldots,\delta_m),
\]

the exact dual basis is

\[
\begin{array}{c|c|c}
\text{wall}&\text{dual vector}&\text{square}\\ \hline
D&c&0\\
\delta_0&d&0\\
\delta_1&c+d-e_1&1\\
\delta_j,\ 2\le j\le m
 &c+2d-e_1-\cdots-e_j&4-j.
\end{array}                                             \tag{2.1}
\]

The independent checker verifies every Kronecker pairing
\(v_i.r_j=\delta_{ij}\); it does not invert or import the producer's Gram
matrix.  These vectors lie in the same future component, as is also clear
from their positive pairing with \(c+d\).

A Lorentzian Coxeter simplex has finite volume exactly when every
deleted-wall vertex is finite or ideal, equivalently when every dual square
in (2.1) is nonnegative.  Hence

\[
                       \boxed{m\le4}.                   \tag{2.2}
\]

The certificate obtains precisely \(m=2,3,4\) over the hostile range
\(2\le m\le48\).  This range is only a regression: formula (2.1) proves the
unbounded statement.

At \(m=5\),

\[
                  v_5=c+2d-e_1-\cdots-e_5,\qquad v_5^2=-1. \tag{2.3}
\]

If all additional wall normals are nonnegative on \(v_5\), then \(v_5\)
remains in the homogeneous chamber cone.  A finite-volume finite-sided
hyperbolic polyhedron has its homogeneous cone generated by future timelike
and future null vertex rays; the convex future cone contains no spacelike
vector.  Thus some new wall \(C\) must have \(C.v_5<0\).  Writing

\[
                         C=ac+bd-\sum\mu_i e_i
\]

gives exactly

\[
             C.v_5=2a+b-\sum_{i=1}^5\mu_i<0.            \tag{2.4}
\]

This is a necessary hidden-incidence condition.  It neither produces the
curve nor proves that a proposed class is a nef facet.

## 3. The full-NS ruling swap defeats strong maximality

Put \(r=c-d\).  Then

\[
 r^2=-2,\qquad s_r(x)=x+(x.r)r,\qquad
 s_r(c)=d,\quad s_r(d)=c,\quad s_r(e_i)=e_i.            \tag{3.1}
\]

Thus \(s_r\) is an integral involutive isometry of \(L\), and it preserves
the chosen positive component because it fixes \(c+d\).  Multiplying the
form by two changes neither the reflection nor its integrality.

The mirror genuinely meets the ample-cone interior.  Start with \(c+d\)
on \(\mathbf P^1\times\mathbf P^1\).  At each blowup, the standard
small-exceptional-coefficient lemma says

\[
                     f^*A-\epsilon E
\]

is ample for all sufficiently small \(\epsilon>0\).  Iteration produces

\[
                A=\pi^*(c+d)-\sum_i\lambda_i e_i,
                \qquad\lambda_i>0,
\]

and \(A.r=0\), because every \(e_i\) is orthogonal to \(r\).  Openness of
the ample cone also shows that the mirror crosses, rather than merely
supports, its interior.

Now assume the complete nef cone is a finite-volume Coxeter chamber \(P\),
with facet group \(W(P)\).  A nonidentity reflection belonging to a Coxeter
group cannot have its mirror cross the interior of its fundamental
chamber.  Therefore

\[
                              s_r\notin W(P).            \tag{3.2}
\]

Nevertheless

\[
                         \langle W(P),s_r\rangle
                         \subset O^+(L)                  \tag{3.3}
\]

is discrete, because \(O^+(L)\) is discrete, and it is generated by
reflections.  It is proper by (3.2).  If maximality is restricted to
finite-covolume or arithmetic reflection groups, the argument still works:
a discrete overgroup of the finite-covolume lattice \(W(P)\) has finite
index, hence finite covolume and the same arithmetic commensurability
class.

This proves the advertised obstruction for every \(m\), before any hidden
curve analysis.  It uses the **full** Neron--Severi host and the literal
complete nef cone.

## 4. Truncator arithmetic, genus, and what is not automatic

The first spacelike dual vector is also an integral root:

\[
 T_5=c+2d-e_1-\cdots-e_5,\qquad
 T_5^2=-1,\quad K.T_5=-1,\quad p_a(T_5)=0.              \tag{4.1}
\]

Its half-space has the correct sign to cut \(v_5\).  Since
\(h^0(\mathcal O(c+2d))=6\), following a prescribed length-five cluster
imposes at most five linear conditions.  Hence \(T_5\) is always an
effective divisor class at \(m=5\), and its total pullback remains effective
after later blowups.  This does not prove irreducibility, smoothness, or
facet status.  The hostile archive now retains \(T_5\) for \(m>5\), which
the original producer archive omitted.

For a primitive reflecting vector in the unimodular lattice \(L\), the
divisibility is one; reflection integrality therefore forces its square,
after primitive normalization, to be \(-1\) or \(-2\).  Write a compatible
horizontal vector as

\[
                    x=\alpha c+\beta d-\sum_i x_i e_i.
\]

Nonnegative pairing with the affine walls and \(D\) is exactly

\[
x_1\ge\cdots\ge x_m\ge0,\qquad
\beta\ge x_1+x_2,\qquad \alpha\ge x_1.                  \tag{4.2}
\]

There is no negative vector with \(\beta=1\).  At the smallest positive
normalized height \(\beta^2/(-x^2)=2\), one must therefore have
\(\beta=2\) and \(x^2=-2\).  Condition (4.2) either permits one coefficient
equal to \(2\), which makes the square nonnegative, or makes every nonzero
\(x_i\) equal to \(1\).  The square equation then says that their number is
\(4\alpha+2\).  This independently recovers

\[
 R_a=ac+2d-e_1-\cdots-e_{4a+2}.                         \tag{4.3}
\]

Direct intersection gives

\[
\begin{aligned}
R_a^2&=-2,&K.R_a&=2a-2,&p_a(R_a)&=a-1,\\
R_a.D&=a-1,&
R_a.R_b&=2(b-a-1)\quad(a<b),\\
R_a.\delta_{4a+2}&=1,
\end{aligned}                                           \tag{4.4}
\]

with all other affine-wall pairings zero.  The hostile checker verifies
these formulas for \(1\le a\le12\) and all pairings among the first nine
roots.  Formula (4.4), rather than that finite range, is the proof.

The jet-row cost is correct as a necessary condition.  Downstairs,

\[
 h^0\bigl(\mathcal O(ac+2d)\bigr)=3a+3.                 \tag{4.5}
\]

The \(4a+2\) cluster conditions can have a nonzero kernel only if at least
\(a\) rows are dependent.  For an irreducible representative of the exact
class \(R_a\), every virtual multiplicity is one.  After the first two
independent jet conditions, the projective system has dimension \(3a\),
while the free cluster has dimension \(4a\); such a curve determines only
finitely many continuing branch directions.  Its incidence image therefore
has codimension at least \(a\).

That last codimension statement does not cover arbitrary reducible
effective divisors.  Excess multiplicity can leave exceptional fixed
components and allow later cluster directions to move.  The repaired
producer retains the universally valid row-dependency condition and states
the codimension bound only for the irreducible exact-multiplicity locus
relevant to a physical facet.

Only \(R_1\) has arithmetic genus zero.  For \(a\ge2\), \(p_a(R_a)>0\), so
those classes cannot be smooth rational walls.  The repaired producer
correctly phrases this as an obstruction to realizing the **whole
first-height family** in the project's physical wall alphabet.  It does
not claim that every \(R_a\) is required as a facet.

## 5. Full host versus restricted and adapted hosts

The theorem cannot be exported silently to a proper sublattice which omits
\(c-d\), or to a proper invariant slice.  In those settings the obstructing
reflection may simply be unavailable.

For an integral half-boundary overhost

\[
                         M=\langle L,B_i/2\rangle,
                         \qquad q_M=2q_L,
\]

one has

\[
 s_{c-d}(B_i/2)
   =B_i/2+{B_i.(c-d)\over2}(c-d).                       \tag{5.1}
\]

Thus every half-boundary must pass the even-pairing test.  An apparent
escape through \((c-d)/2\in M\) does not weaken this statement: if both
\((c-d)/2\) and \(B_i/2\) lie in the integral doubled host, their pairing
\(B_i.(c-d)/2\) is already integral.  Hence odd pairing is impossible in
that case as well.

The full-NS obstruction therefore remains exact, while an adapted-host
extension requires its own parity calculation.  The producer and this
audit correctly make no universal adapted-host claim.

## 6. Scope of this audit

The companion affine-H3 abundance discussion in Section 5 of the producer
is not used anywhere in the affine-\(B_m\) obstruction.  Its semiampleness
and index-one claims have a separate independent audit in the
[affine-H3 hostile audit](rational_affine_h3_finite_volume_abundance_hostile_audit.md).
The present audit deliberately keeps the two arguments logically
independent.

The final boundary is therefore sharp:

- the explicit affine packet is not an unbounded full-NS
  reflection-maximal counterfamily;
- its universally visible simplex stops being finite-volume at \(m=5\);
- any nonmaximal or relative-host completion from \(m=5\) onward owes
  hidden horizontal incidence;
- no claim excludes specially tuned finite-volume nonmaximal nef cones;
- no root census or external classification is being substituted for the
  algebraic-geometric argument.

The source artifacts pinned by the hostile certificate are:

    0c3e32ef6763308cfbe4b4b28e21452c428611a9956d0a5a9333220e55d69b60  results/rational_affine_b_tail_global_completion_obstruction.md
    a76414589ce0c7680fcee6ea1928536e540d324189b456eb834cc125c2efd517  scripts/certify_rational_affine_b_tail_global_completion_obstruction.py
    996d0919993e0b9756f1cbf36c6e922a2e4fa084cf4ad86643326253768fbc38  data/enumerations/rational_affine_b_tail_global_completion_obstruction.json
    5582bfc90b256d5f5dc91aa6622c3e9ee1e2cc200f1e77b59fe218acac36f500  tests/test_rational_affine_b_tail_global_completion_obstruction.py
