# Physical minus-two narrow stars force finitely many compact K3 trace hosts

**Status (2026-08-11): proved, conditionally on an exact physical-star
lift.**  This is a strengthening of the
[compact narrow-facet field theorem](k3_compact_narrow_facet_field_finiteness.md).
It does not prove that the required physical lift exists.

Let \(P\subset\mathbf H^{m-1}\), \(m\ge3\), have an exact compact
primitive full-trace complex-K3 carrier.  Choose the minimal-aperture facet
\(E\) used in the field theorem.  Add the following genuinely geometric
hypothesis:

> Every normal in the star
> \[
>  \mathcal S(E)=\{E\}\cup\{F:F\cap E\ne\varnothing\}
> \]
> is represented by the **same vector** in the full trace module as a
> primitive physical smooth-rational curve class of square \(-2\).

“Same vector” is essential.  It is not enough that an unrelated divisor
have the same projected wall line in the distinguished real eigenspace.

Under this hypothesis there are only finitely many marked integral
full-trace hosts, including all integral, even, order-stable saturations and
overlattices.  The proof is classification-free.  Its new mechanism is the
placewise energy identity

\[
 \sum_{\tau:k\hookrightarrow\mathbf R}
       -\tau(g(\delta_i,\delta_i))=2,                 \tag{0.1}
\]

where \(g\) is the \(k\)-valued trace kernel of the K3 intersection form.
The minimal-aperture bound then puts every \(g(\delta_i,\delta_j)\) in a
bounded box of an order codifferent.  Integral Cartan cycles bound the
coefficient order from below, and ordinary trace duality bounds every
possible module from above.

If, in addition, reflection maximality is taken relative to the **same
marked \(k\)-reflection host**, then only finitely many compact maximal
chambers occur.  For \(k\ne\mathbf Q\), this last statement cannot be
obtained by silently replacing \(k\)-reflections with ordinary rank-one
K3 root reflections.

## 1. Exact lattice package

Let \(k\) be totally real with distinguished embedding \(\tau _0\), let
\(\mathcal O\subset\mathcal O_k\) be an order, and let \(M\) be a finite
projective rank-\(m\) \(\mathcal O\)-module.  Write

\[
                         V=M\otimes_{\mathcal O}k.
\]

Let \(N=M_{\mathbf Z}\) carry an even integral form \(B\), and suppose the
\(\mathcal O\)-action is self-adjoint.  Equivalently, there is a unique
symmetric \(k\)-bilinear trace kernel \(g\) such that

\[
 B(ax,y)=\operatorname {Tr}_{k/\mathbf Q}
                 \bigl(a\,g(x,y)\bigr)
 \quad (a\in k,\ x,y\in V).                            \tag{1.1}
\]

The sign convention is the K3 convention:

\[
 \operatorname {sig}(g^{\tau _0})=(1,m-1),
 \qquad g^\tau<0\quad(\tau\ne\tau _0).                 \tag{1.2}
\]

Thus \(B\) has signature \((1,dm-1)\), where
\(d=[k:\mathbf Q]\).  Assume \(N\) embeds primitively in the Neron--Severi
lattice of a complex projective K3 surface.  Hence

\[
                              dm\le20.                 \tag{1.3}
\]

The facet reflection attached to a normal \(r\in V\) is

\[
 s_r(x)=x-\frac{2g(x,r)}{g(r,r)}\,r.                   \tag{1.4}
\]

The exact carrier hypothesis says that all displayed facet reflections
preserve \(M\).  The chamber in the distinguished eigenspace is compact and
is exactly the supplied \(P\).

Choose a minimal-aperture facet \(E\).  The compact theorem proves two facts:

1. all absolute normalized pairings in \(\mathcal S(E)\) are strictly smaller
   than \(14\);
2. the star spans \(V\), so it contains \(m\) normals forming a \(k\)-basis.

The extra hypothesis in this note lets us choose that basis as actual trace
vectors

\[
                  \delta _1,\ldots,\delta _m\in M,
 \qquad B(\delta_i,\delta_i)=-2.                       \tag{1.5}
\]

Each \(\delta_i\) is the class of the physical curve defining that very
wall.  Primitivity is automatic from square \(-2\) in an integral lattice,
but is retained to make the geometry explicit.

## 2. The physical square is a placewise energy budget

Put

\[
 g_{ij}=g(\delta_i,\delta_j),\qquad
 t_{i,\tau}=-\tau(g_{ii}).
\]

At \(\tau_0\), a facet normal is negative, and at every other real place
the whole form is negative definite.  Therefore

\[
                         t_{i,\tau}>0.                 \tag{2.1}
\]

Taking the trace in (1.5) gives the exact identity

\[
                         \sum_\tau t_{i,\tau}=2.       \tag{2.2}
\]

In particular \(t_{i,\tau}\le2\) at every place.

For \(i\ne j\), the narrow-star estimate at the distinguished place says

\[
 \frac{2\lvert\tau_0(g_{ij})\rvert}
 {\sqrt{t_{i,\tau_0}t_{j,\tau_0}}}<14,
\]

and hence

\[
 \lvert\tau_0(g_{ij})\rvert<14.                       \tag{2.3}
\]

At a nonidentity place, negative definiteness and Cauchy--Schwarz give

\[
 \lvert\tau(g_{ij})\rvert
 \le\sqrt{t_{i,\tau}t_{j,\tau}}\le2.                  \tag{2.4}
\]

Thus every conjugate of every entry of the \(k\)-Gram matrix
\(G=(g_{ij})\) lies in a fixed compact interval.  One also gets the useful
integral-trace check

\[
\begin{aligned}
 \lvert B(\delta_i,\delta_j)\rvert
 &\le
 7\sqrt{t_{i,\tau_0}t_{j,\tau_0}}
 +\sum_{\tau\ne\tau_0}\sqrt{t_{i,\tau}t_{j,\tau}}\\
 &<14+2=16.                                           \tag{2.5}
\end{aligned}
\]

The second term is at most \(2\) by Cauchy--Schwarz applied to the vectors
of the remaining \(t\)-coordinates.  Equation (2.5) is not by itself
enough to recover the field form; the codifferent step below is what keeps
all \(k\)-coordinates.

## 3. Integral cycles leave only finitely many coefficient orders

For the selected star basis put

\[
                         a_{ij}=\frac{2g_{ij}}{g_{ii}}
\]

For \(T_i=1-s_{\delta_i}\), rank-one calculus on the finite projective
\(\mathcal O\)-module gives

\[
 \operatorname {tr}_{\mathcal O}
 (T_{i_1}\cdots T_{i_\ell})
 =a_{i_1i_2}\cdots a_{i_\ell i_1}\in\mathcal O.       \tag{3.1}
\]

The compact field theorem supplies a connected Lorentz block whose cyclic
Cartan products generate all of \(k\).  Closed-walk and chord splitting
leave a finite packet \(\mathcal C\) of pair cycles and simple chordless
cycles, each on at most \(m\) nodes.  Every \(c\in\mathcal C\) is an algebraic
integer and

\[
                         \lvert\tau(c)\rvert<14^m      \tag{3.2}
\]

at every real embedding.

There are only finitely many possible packets \(\mathcal C\).  Indeed,
(1.3) bounds \(m\le20\); if \(d>1\), it improves this to
\(d\le6\) and \(m\le10\), while \(d=1\) gives ordinary bounded integers.
Bounded degree and bounded conjugates put the coefficients of every monic
minimal polynomial in finite integer boxes.

For each packet set

\[
                         \mathcal O_0=\mathbf Z[\mathcal C].  \tag{3.3}
\]

Because the fraction field generated by \(\mathcal C\) is \(k\),
\(\mathcal O_0\) is a full order in \(k\).  Equation (3.1) gives

\[
                 \mathcal O_0\subseteq\mathcal O
                              \subseteq\mathcal O_k.   \tag{3.4}
\]

The quotient \(\mathcal O_k/\mathcal O_0\) is finite.  Every intermediate
order is, additively, the inverse image of a subgroup of that finite
quotient; imposing multiplication closure only removes candidates.
Consequently only finitely many coefficient orders \(\mathcal O\) occur.

This is the point that forbids arbitrary conductors.  Merely knowing that
\(k\) belongs to a finite field set would not do so.

## 4. The codifferent box leaves finitely many field Gram matrices

For an order \(\mathcal O\), write

\[
 \mathcal O^\vee=
 \{z\in k:\operatorname {Tr}_{k/\mathbf Q}(z\mathcal O)
                      \subseteq\mathbf Z\}.            \tag{4.1}
\]

This is the order codifferent, a full lattice in the Minkowski space
\(k\otimes_{\mathbf Q}\mathbf R\).

Since \(a\delta_i,\delta_j\in M\) for every \(a\in\mathcal O\), integrality
of \(B\) and (1.1) give

\[
 \operatorname {Tr}_{k/\mathbf Q}(a g_{ij})
       =B(a\delta_i,\delta_j)\in\mathbf Z.
\]

Therefore

\[
                              g_{ij}\in\mathcal O^\vee. \tag{4.2}
\]

For each of the finitely many orders, (2.2)--(2.4) put every diagonal and
off-diagonal entry in the intersection of \(\mathcal O^\vee\) with a fixed
compact Minkowski box.  A lattice has finite intersection with a compact
set.  Hence only finitely many symmetric matrices \(G=(g_{ij})\) occur.

This absorbs all choices of a conventional trace scalar.  Relative to the
sign convention of the compact field theorem one writes

\[
                              g=-\alpha h
\]

with \(\alpha\) totally positive.  Presentations can vary by the simultaneous
rescaling \((\alpha,h)\mapsto(\alpha/u,uh)\), but the trace host and the
native reflection

\[
 x\longmapsto x-\frac{2g(x,r)}{g(r,r)}\,r
             =x-\frac{2h(x,r)}{h(r,r)}\,r
\]

depend only on \(g\), not on \(\alpha\) and \(h\) separately.  The matrix
\(g\) is now in a finite set.

## 5. Root-frame duality leaves finitely many modules and saturations

Fix one surviving triple \((k,\mathcal O,G)\), and let

\[
                     F=\bigoplus_{i=1}^m\mathcal O\delta_i
                         \subseteq M.                  \tag{5.1}
\]

The \(\delta_i\) are \(k\)-independent, so \(F\) is a free
\(\mathcal O\)-module of rank \(m\) and has finite additive index in \(M\).
Choose a \(\mathbf Z\)-basis
\(\omega_1,\ldots,\omega_d\) of \(\mathcal O\).  The complete integral
Gram matrix of \(F\) is

\[
 C_{(i,a),(j,b)}
   =\operatorname {Tr}_{k/\mathbf Q}
          (\omega_a\omega_b g_{ij}).                   \tag{5.2}
\]

Thus \(F\), not merely its normalized Coxeter diagram, belongs to a finite
effective list.

Let

\[
 F^\#=\{x\in F\otimes\mathbf Q:B(x,F)\subseteq\mathbf Z\}.
\]

Because \(F\subseteq M\) and \(M\) is integral,

\[
                              F\subseteq M\subseteq F^\#.     \tag{5.3}
\]

The finite discriminant group satisfies

\[
 |F^\#/F|=|\det C|.                                    \tag{5.4}
\]

Hence \(M/F\) is one of finitely many subgroups of \(F^\#/F\).  More
quantitatively,

\[
 [M:F]^2\mid |\det C|                                  \tag{5.5}
\]

because
\(\det F=[M:F]^2\det M\) and \(\det M\) is a nonzero integer.
Order stability, projectivity, evenness, preservation by the displayed
reflections, and primitive K3 embeddability are all finite filters on these
intermediate groups.

The same argument includes any integral saturation \(N'\) of the
root-generated module in its rational span:

\[
                              F\subseteq N'\subseteq F^\#.
\]

For even overlattices, \(N'/F\) is an isotropic subgroup of the finite
discriminant quadratic module of \(F\); requiring \(\mathcal O\)- and
reflection-stability again only removes candidates.  No ideal-class or
overlattice classification is needed for finiteness.

### Theorem 5.1 (finite physical-star trace hosts)

Across all exact compact primitive full-trace complex-K3 carriers in
hyperbolic dimension at least two satisfying the physical-star hypothesis
of Section 1, only finitely many isomorphism classes of marked
\((k,\tau_0,\mathcal O,M,g)\) trace hosts and their integral even
saturations occur.  After adjoining a choice \(\mathcal C^+\) of component of
the positive cone of \(g^{\tau_0}\), only finitely many complete oriented
packets

\[
 \mathcal H=(k,\tau_0,\mathcal O,M,g,\mathcal C^+)     \tag{5.6}
\]

occur as well: a field has finitely many real embeddings and a Lorentz
positive cone has two components.  The list is effective from the finite
cycle, overorder, codifferent-box, and discriminant-subgroup enumerations
above.

### Corollary 5.2 (exact full-Picard eigen-nef return)

Suppose the exact carrier has

\[
                         \operatorname {NS}(X)=N
\]

and returns \(P\) as the exact complete eigen-nef intersection

\[
             P=\mathbf P\bigl(\operatorname {Nef}(X)\cap E_{\tau_0}\bigr)
\]

with every eigenslice facet exposed by an actual K3 nef wall.  Such a wall
is generated by a primitive smooth-rational square-\(-2\) curve class
\(\delta\).  Since \(\operatorname {NS}(X)=N=M_{\mathbf Z}\), one has
\(\delta\in M\), and its distinguished projection is the native wall
normal.

Here is the elementary descent of the proportionality scalar.  In a
\(k\)-basis write the arithmetic normal as \(r=(r_j)\) and the curve vector
as \(\delta=(\delta_j)\).  Equality of their distinguished real lines gives
\(\tau_0(r)=c\,\tau_0(\delta)\).  Choose \(\delta_\ell\ne0\) and put
\(a=r_\ell/\delta_\ell\in k\).  Then \(c=\tau_0(a)\), and injectivity of
\(\tau_0:k\hookrightarrow\mathbf R\) applied coordinatewise gives
\(r=a\delta\).  Thus rescaling the former to the primitive curve vector
does not change its native reflection.  The physical-star hypothesis is
therefore automatic in this full-Picard exact eigen-nef lane, and
Theorem 5.1 applies.

For \(k=\mathbf Q\), the distinguished eigenspace is the whole trace
space, so this specializes to literal full-nef return.  For
\(k\ne\mathbf Q\), one must retain the eigen-nef formulation: literal
full-nef equality would have the wrong dimension and would force the
rational lane.

The equality \(\operatorname {NS}(X)=N\) is doing real work here.  If the
trace block is merely primitive in a larger Picard lattice, a physical wall
class can have a nonzero component outside the trace block; equality of its
projected wall with a field-normal line does not put that curve vector in
\(M\).

## 6. Safe composition with same-host maximality

For the complete oriented host (5.6), define

\[
 \Gamma_{\mathcal H}=\{u\in\operatorname {GL}_{\mathcal O}(M):
   g(ux,uy)=g(x,y),\
   u_{\tau_0}(\mathcal C^+)=\mathcal C^+\}
\]

and let \(R(\mathcal H)\) be the subgroup generated by every element of
\(\Gamma_{\mathcal H}\) whose distinguished action is a hyperplane reflection.
This is exactly the complete native host and reflection species of the
independently frozen
[compact trace-host composition audit](compact_trace_host_chamber_finiteness_composition_hostile_audit.md).
Admissibility makes the other real orthogonal factors compact, so the
distinguished image of \(\Gamma_{\mathcal H}\) is discrete and its complete
native mirror arrangement is locally finite.

Assume the displayed facet group \(W(P)\) is reflection-maximal relative to
this same marked host.  If an integral \(k\)-reflection lay outside
\(W(P)\), adjoining it would give a strictly larger discrete group generated
by reflections inside \(\Gamma_{\mathcal H}\).  Theorem 2.1 of the composition
audit therefore gives

\[
                              W(P)=R(\mathcal H).        \tag{6.1}
\]

That theorem also proves that the locally finite complete arrangement has
one \(R(\mathcal H)\)-orbit of chambers.  Since Theorem 5.1 gives finitely many
complete oriented packets, it consequently gives:

### Corollary 6.1 (finite same-host maximal chambers)

There are only finitely many isometry classes of compact chambers satisfying
Theorem 5.1 and reflection maximality relative to their marked integral
\(k\)-reflection host.

This formulation is necessary.  If \(d>1\), a \(k\)-reflection has a
\(d\)-dimensional minus eigenspace after restriction of scalars to
\(\mathbf Q\).  The ordinary K3 root reflection in the physical class
\(\delta_i\) has a one-dimensional rational minus eigenspace.  They are not
the same operator, even though the distinguished eigenspace sees the same
facet wall.  Thus maximality for ordinary rank-one reflections of the full
Neron--Severi lattice cannot replace the hypothesis of Corollary 6.1.

## 7. Exact scope and the remaining geometric hypothesis

The conclusion is about the full trace/eigen reflection host and every
integral saturation in its rational span.  It does not assert finiteness of
arbitrary ambient Neron--Severi extensions perpendicular to the trace block,
nor finiteness of carrier surfaces or periods.  Full-Picard equality or a
separate bounded-complement theorem would be needed for that stronger
claim.

The physical-star premise is precisely the nontrivial geometric information
missing from a bare \(\operatorname {NS}\)-embedding:

1. the star must lift by the same trace-module vectors, not merely by
   proportional distinguished projections;
2. those vectors must be simultaneously effective irreducible
   square-\(-2\) curves;
3. the selected physical star must span the \(k\)-module.

Without the square-\(-2\) identity, multiplying a trace form by an
unbounded positive integer leaves every normalized wall pairing and every
field reflection unchanged while its integral discriminant escapes.
Without a physical \(k\)-basis, uncontrolled directions of the trace form
remain.  The theorem therefore does not prove Property 4, but it shows that
this explicit physical-star property is sufficient to upgrade the frozen
field theorem to finite trace hosts and, with same-host \(k\)-reflection
maximality, to a finite compact chamber set.
