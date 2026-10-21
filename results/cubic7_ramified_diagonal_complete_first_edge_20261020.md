# Complete first edge of the ramified-seven diagonal family

This paper proof uses no new computation. Let
\(t=2\cos(2\pi/7)\), \(t^3+t^2-2t-1=0\),
\(\beta=2-t\), \(\alpha=t\beta\),
\(u=t^2+t-1\), and \(T_n=\langle-\alpha\rangle\oplus I_n\)
over \(O=\mathbf Z[t]\), with \(n\ge5\).
The old vertex is \(p=(1;0)\), with outward simple roots
\(\alpha_i=-e_i+e_{i+1}\) for \(i<n\), \(\alpha_n=-e_n\).
We release \(\alpha_1\) and retain \(B_{n-1}\).

The complete primitive alphabet \(\{1,2\}\), full \(B_n\) zero germ,
coordinate-facet identity, and exact I2(7) section/glue interface are imported
from [the separately reviewed reduction](cubic7_ramified_diagonal_norm_and_facet_reduction_20261020.md),
SHA256 `34cef661bd0333b633d4cf0912b096a273635d3bff54db8fa358b7222385e87a`.
Its field and trace inputs refer to
[the arithmetic calibration](cubic7_diagonal_first_edge_and_literature_gate_20261020.md),
SHA256 `f6251aacd1259f43527c6aa6ebd37c579edfbb0f4b5a1349c517545f9d6ac1f6`.
No norm/face proof from those sources is repeated here. Their literature
boundary is retained: no exact published theorem excluding this
nonunimodular family has been established by the bounded screen.

## Complete first finite edge, in every dimension at least five

Let \(a=t-1\), a unit. The identity \(\alpha=1-a^2\) gives
\(X=a^{-1}\in O\). For \(n\ge2\),

\[
v=(X;X,0,\ldots),\qquad d=(X;X,1,0,\ldots),\qquad
P=(X;\alpha X,0,\ldots)
\]

have norms one, two, and \(-\alpha\), respectively. Both mirrors meet the
edge retaining \(B_{n-1}\) at \(P\); \(d\) is acute to all its retained
walls. In the binary coordinates, the integral involution

\[
J=\begin{pmatrix}X&-X\\\alpha X&-X\end{pmatrix}
\]

preserves the Gram matrix, sends \(p\) to \(P\), and sends the old
\(B_n\) simple root \(-e_1+e_2\) to \(d\), fixing the other simple roots.
This is a finite upper, but it will not be the first endpoint. We now
exhaust all possible earlier finite continuations.

Retain \(B_{n-1}\) on coordinates \(2,\ldots,n\). A finite rank-one
continuation is either an orthogonal root of norm one or two, or the
long-end \(B_n\) continuation of norm two. (A finite system containing
\(B_{n-1}\), of rank at least four, permits no other connected extension.)
The latter has spatial coordinates \((y,1,0,\ldots)\); the former has
\((y,0,\ldots)\). Thus the scalar equation is

\[
y^2-\alpha x^2=q_0,\qquad q_0=1\text{ or }2,
\]

where the connected case has residual \(q_0=1\). Choose the inward time
orientation \(x>0\), and put \(A=ax\). A collision no later than the upper
has \(y\ge x\), so \(A^2\le q_0\). At the other two embeddings,

\[
A_i^2\le q_0\frac{a_i^2}{a_i^2-1}.
\]

The polynomial of \(\alpha\) is \(z^3+7z^2-7\), so
\(\operatorname{Tr}(\alpha^{-1})=0\). It follows that

\[
\operatorname{Tr}(A^2)\le q_0(3+\alpha^{-1})
 <\frac{61q_0}{15},
\]

because \(1<t<5/4\) gives \(\alpha>15/16\).
For \(A=A_0+A_1t+A_2t^2\), the exact trace form is

\[
\operatorname{Tr}(A^2)=
\frac{(3A_0-A_1+5A_2)^2+
14(A_1^2-A_1A_2+A_2^2)}3.
\]

For \(q_0=1\), trace at most four forces \(A=1\). This gives
\(x=y=X\), hence only the displayed upper roots \(v,d\).

For \(q_0=2\), trace at most eight gives the signed list

\[
0,\ 1,\ t,\ v_0=t^2-2,\ u,\ t+1,\ v_0+1,\ u-1.
\]

Indeed the hexagonal form \(A_1^2-A_1A_2+A_2^2\) is zero or one; in
the latter case the remaining square is at most ten. The physical
condition \(0<A\le\sqrt2\) leaves
\(1,t,-v_0,v_0+1,u-1\).

For complete rational conjugate filters, order the other embeddings as
\(t_2=v_0\in(-1/2,-2/5)\), \(t_3=-u\in(-2,-7/4)\).
Their necessary bounds are
\(A_2^2<49/12\) and \(A_3^2<242/105\).
The third conjugate of \(t\) has square \(u^2=t+2>3\), excluding it;
the third conjugate of \(v_0+1\) is \(t+1>11/5\), excluding it; and the
second conjugate of \(u-1\) is \(-t-1<-11/5\), excluding it.
Only \(A=1,-v_0\) remain.

For \(A=1\), the equation requires \(y^2=X^2+1\). This cannot be a field
square since

\[
N(X^2+1)=N((t-1)^2+1)=N(t^2-2t+2)=29.
\]

For \(A=-v_0=2-t^2\), the identities
\((2-t^2)X=u\) and \(\alpha u^2+2=(t+1)^2\) give exactly

\[
\boxed{r=(u;t+1,0,\ldots),\quad r^2=2.}
\]

Its collision is strictly earlier than the upper, since \(u<t+1\).
The endpoint is

\[
Q=(t+1;\alpha u,0,\ldots),\qquad Q^2=-2\alpha.
\]

There is no simultaneous second simple wall: the above inventory exhausts
the finite extensions, and the only other candidate collisions are later.
Local finiteness of the integral mirror arrangement on the compact edge
segment now proves that \(r\) is the actual first wall. The full new germ
is \(B_{n-1}+A_1\); an enlargement to \(B_n\) would have supplied the
excluded earlier connected continuation.

Finally \(u(t+1)=X\) and \((t+1)^2-1=X\) show that the upper involution
\(J\) is exactly \(s_r\). It is an actual reflection, NOT a chamber
stabilizer or an infinite-period certificate. Its formal endpoint \(P\)
is the reflected old vertex beyond \(r\). No reflectivity decision follows
from this first edge alone.

The literal rank-four member already has bounded archived prefixes,
[without a nonreflectivity conclusion](cubic7_timelike_scaled_rank5_exact_prefixes.md).
The complete first edge above is an own paper argument, not a claim that
that bounded search terminates. No new example, full-host exclusion,
worker or scientific execution is claimed.
