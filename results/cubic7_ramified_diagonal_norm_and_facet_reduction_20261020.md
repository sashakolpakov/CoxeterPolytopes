# The ramified-seven diagonal family: exact norm and face reductions

This is a paper argument, not a new Vinberg computation or a reflectivity
claim. Put

\[
t=2\cos(2\pi/7),\quad t^3+t^2-2t-1=0,\quad
\beta=2-t,\quad \alpha=t\beta,
\qquad T_n=\langle-\alpha\rangle\oplus I_n
\]

over \(O=\mathbf Z[t]\), at the embedding \(1<t<5/4\).
The other two conjugates of \(t\) are negative, whereas all three conjugates
of \(\beta\) are positive. Thus this is an admissible hyperbolic lattice.
Its determinant ideal is \((\beta)\), of norm seven; it is not integrally
isometric to the unimodular diagonal family.

The field and trace facts used below are proved in
[the first-edge calibration](cubic7_diagonal_first_edge_and_literature_gate_20261020.md):
\(O\) is the full ring, its class number is one, all totally positive units
are squares, and 2 is inert. For nonzero \(z\in O\),
\(\operatorname{Tr}(z^2)\ge3\), with equality only for \(z=\pm1\).

## The apparent four-class norm pool reduces to two classes

For a primitive root \(r=(x;z_1,\ldots,z_n)\), write \(q=r^2\).
Reflection integrality is exactly

\[
q\mid2(-\alpha x,z_1,\ldots,z_n).
\]

Outside 2 and \(\beta\), primitivity makes the pairing ideal a unit.
At 2 the lattice is unimodular locally, so \(v_2(q)\le1\).
At \(\beta\), either a spatial coordinate is a unit and \(v_\beta(q)=0\),
or all spatial coordinates are divisible by \(\beta\), the time coordinate
is a unit, and
\(q\equiv-t\beta x^2\pmod{\beta^2}\) gives \(v_\beta(q)=1\).
Unit normalization therefore gives the conservative list
\(1,2,\beta,2\beta\).

**Neither ramified class occurs, in any rank.** If \(q=\epsilon\beta\),
\(\epsilon\in\{1,2\}\), reflection integrality gives \(z_i=\beta w_i\).
The norm equation becomes

\[
\beta\sum_iw_i^2-tx^2=\epsilon.
\]

But \(O/(\beta)=\mathbf F_7\) and \(t=2\) there. Thus \(x^2=3\) for
\(\epsilon=1\), or \(x^2=6\) for \(\epsilon=2\). Both are nonsquares
in \(\mathbf F_7\). Consequently the complete primitive norm alphabet is
\(\boxed{\{1,2\}}\), although the lattice is nonunimodular.

At \(p=(1;0)\), a zero-height root has only spatial coordinates. Norm one
gives precisely \(\pm e_i\); norm two gives precisely \(\pm e_i\pm e_j\).
Indeed at most two coordinates can be nonzero by the trace bound; equality
forces each to be \(\pm1\), while a single coordinate of square two would
put \(\sqrt2\) in a cubic field. The FULL zero-height germ is therefore
the ordinary \(B_n\), not merely a selected subsystem.

## Full coordinate-facet identity

The coordinate wall \(e_n^\perp\) is an actual wall at this germ. A norm-one
mirror crossing it nonorthogonally would have last coordinate \(b\in O\)
with \(b_i^2\le1\) at all embeddings and strict inequality at the
hyperbolic embedding, contrary to \(|N(b)|\ge1\).

For a norm-two crossing mirror, the last coordinate satisfies
\(b_i^2\le2\), with one strict. Thus \(\operatorname{Tr}(b^2)\le5\).
The exact trace form in the calibration shows that the possibilities are
\(0,\pm1\) and the signed conjugates of \(t\); each conjugate orbit of
\(t\) has a member of magnitude \(2\cos(\pi/7)>\sqrt2\), so only
\(b=0,\pm1\) remain. Orthogonal projection subtracts \(b e_n\) and gives
an integral root of norm two or one in \(T_{n-1}\).

Conversely every primitive norm-one or norm-two root of \(T_{n-1}\)
extends by a zero coordinate to a reflector of \(T_n\). Therefore the
intrinsic arrangement on this actual facet is EXACTLY the full integral
root arrangement of \(T_{n-1}\). In particular, a compact full chamber in
dimension \(n\) forces a compact full chamber in every smaller dimension.
This statement is a reduction, not an established low-dimensional
nonreflectivity result.

## The actual I2(7) section of the unimodular family

Put \(u=t^2+t-1\), so \(u^2=t+2\). In the unimodular family
\(\langle-t\rangle\oplus I_{n+2}\), the pair

\[
a_1=(0;-1,1,0,\ldots),\qquad r=(1;u,0,\ldots)
\]

has type \(I_2(7)\). Its integral orthogonal complement has basis
\((u;t,t,0,\ldots),e_3,\ldots,e_{n+2}\), and Gram \(T_n\).
The direct sum with the selected pair has determinant-ideal index
\((\beta)\), of abelian cardinality seven, so it is not a unimodular
splitting. Nevertheless every root of \(T_n\) has norm one or two and
extends to an integral ambient reflector in this basis.

Whenever the displayed pair bounds an ACTUAL face of the ambient chamber
(as in the calibrated first edge in ambient dimension at least five),
every actual neighboring wall is perpendicular to both selected walls:
its finite rank-three diagram contains \(I_2(7)\), and the classification
of finite rank-three systems permits only \(I_2(7)+A_1\).
Hence this actual section chamber is the FULL \(T_n\) chamber. The possible
extra \(\beta,2\beta\) mirrors that would otherwise obstruct this claim
were excluded above. Nonreflectivity of the larger ambient chamber alone
does not imply nonreflectivity of this section.

## A finite upper for a future first-edge census

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
This proves a genuine finite upper and a possible full \(B_n\) endpoint.
It does NOT prove that no earlier wall intervenes, nor an infinite period.

## Existing evidence and literature boundary

The literal \(T_4\) already appears in
[the rank-five archive](cubic7_timelike_scaled_rank5_exact_prefixes.md),
with bounded exact prefixes through height1000. That report explicitly
does not claim nonreflectivity. Its four-class candidate list was
conservative; the modular argument above removes the two empty classes.
The E8-plus-one coefficient archive uses a different positive module.

The bounded primary-literature screen has not produced an exact theorem
excluding this nonunimodular family. The Bugaenko1992 unimodular
classification is not such a theorem without an additional argument.
There is no new compact example, reflectivity decision, worker, or
scientific execution in this note.
