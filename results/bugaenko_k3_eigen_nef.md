# Bugaenko's Chamber as an Exact K3 Eigen-Nef Section

Status: exact existential result for a marked K3 surface.  The arithmetic and
wall enumerations are executable certificates.  The period-map and Weyl-cone
steps use standard K3 theorems and are identified separately below.  No
explicit projective equation or explicit Weyl word is claimed.

## The geometric object

Let

\[
F=\mathbb Q(\sqrt5),\qquad \mathcal O_F=\mathbb Z[\phi],\qquad
M=\langle-2\phi\rangle\oplus E_8,
\]

and let \(h\) be Bugaenko's \(F\)-bilinear form.  On the underlying
rank-18 abelian group put

\[
B(x,y)=\operatorname{Tr}_{F/\mathbb Q}
\left(-\frac{\phi}{\sqrt5}h(x,y)\right).
\]

The earlier trace certificate proves

\[
N=(M,B)\simeq
\langle2\rangle\oplus\langle-2\rangle\oplus E_8(-1)^2
\]

and gives a primitive embedding \(N\hookrightarrow\Lambda_{K3}\).  Its
orthogonal complement has signature \((2,2)\), so the corresponding period
domain has complex dimension two.  Outside the countable union of rational
hyperplane sections defined by classes in \(\Lambda_{K3}\setminus N\), a
period has Neron--Severi lattice exactly \(N\).  Surjectivity of the K3 period
map therefore gives a projective K3 surface \(X\) with
\(\operatorname{NS}(X)\simeq N\).  This is the standard generic-period
argument in lattice-polarized K3 theory; see
[Dolgachev](https://arxiv.org/abs/alg-geom/9502005).

Multiplication on \(M\) supplies a ring embedding

\[
\iota:\mathcal O_F\longrightarrow
\operatorname{End}_{\mathbb Z}(\operatorname{NS}(X)).
\]

It is self-adjoint because

\[
B(ax,y)=\operatorname{Tr}(\alpha a h(x,y))=B(x,ay).
\]

The executable certificate checks both
\(\iota(\phi)^2=\iota(\phi)+1\) and self-adjointness in the integral trace
basis.  It also checks that \(\iota(\phi)\) is **not** an isometry of \(B\).
Thus the object is precisely

\[
(X,\iota:\mathcal O_F\hookrightarrow
\operatorname{End}_{\mathbb Z}\operatorname{NS}(X)),
\]

not an \(\mathcal O_F\)-action by automorphisms of \(X\), and not a claim of
natural real multiplication on the full K3 Hodge structure.

If desired, every rational endomorphism of \(\operatorname{NS}(X)_\mathbb Q\),
including \(\iota(\phi)\), can be represented by a rational algebraic
correspondence on \(X\times X\): products \(D_i\times D_j\) act by
\(x\mapsto(x\cdot D_i)D_j\), and the inverse intersection matrix produces
the dual basis.  This elementary observation gives a class in
\(\operatorname{CH}^2(X\times X)_\mathbb Q\); it does not make the action
canonical, integral as a correspondence, or induced by automorphisms.

## The five K3 shells collapse to norm-two mirrors

Let \(E\subset N_\mathbb R\) be the eigenspace of \(\iota(\phi)\) for the
identity embedding \(\phi\mapsto(1+\sqrt5)/2\).  If
\(\delta\in N=M\) has \(B(\delta,\delta)=-2\) and its hyperplane meets the
positive cone of \(E\), then

\[
h(\delta,\delta)\in
\{5-3\phi,\ 4-2\phi,\ 3-\phi,\ 2,\ 1+\phi\}.
\]

Bugaenko's field Gram matrix is even over \(\mathcal O_F\): its diagonal is
in \(2\mathcal O_F\), while every off-diagonal contribution to
\(h(x,x)\) already occurs twice.  Hence
\(h(x,x)\in2\mathcal O_F\) for every \(x\in M\).  Of the five shells only

\[
4-2\phi\quad\hbox{and}\quad2
\]

belong to \(2\mathcal O_F\).  They define the same hyperplanes because

\[
\phi^2(4-2\phi)=2,
\qquad (\phi-1)^2\,2=4-2\phi,
\]

and both \(\phi\) and \(\phi-1\) are units.  Conversely every integral
norm-two vector is an integral field-reflection root: for
\(h(r,r)=2\), its reflection coefficients are simply the entries of \(Gr\).
Therefore

> the full arrangement obtained by restricting all K3 \((-2)\)-hyperplanes
> to \(E\) is exactly Bugaenko's norm-two reflection arrangement.

Each restricted mirror has exactly four signed trace-root representatives:
\(\pm r\) in the norm-two shell and \(\pm(\phi-1)r\) in the
\(4-2\phi\) shell.  These are two distinct K3 hyperplanes in the ambient
rank-18 positive cone, although their restrictions to \(E\) coincide.
Indeed, equality of two \(F\)-rational hyperplanes in \(E\) makes their
normals proportional by a scalar \(c\in F\); the two norm equations give
\(c^2=1\) or \(c^2=2-\phi=(\phi-1)^2\), respectively.  Thus no further
vectors are hidden behind the same restricted mirror.

## No extra wall cuts the chamber

The audited Vinberg transcript exhausts the simple norm-two roots, and the
independent cone certificate proves compact termination.  Its chamber \(P\)
has 11 facets and 41 finite vertices.  Since the restricted K3 arrangement is
the same norm-two arrangement, no K3 wall cuts the interior of \(P\), and its
only restricted facets are exactly the eleven Bugaenko facets.  In particular
there is no additional non-reflection divisor wall.

There are nevertheless norm-two mirrors which touch lower-dimensional faces
of the **closed** chamber.  They can be enumerated without a norm or radius
cutoff.  A mirror meeting \(\overline P\) cannot meet \(P^\circ\), so its
contact is a nonempty face.  Compactness makes that face contain a vertex.  At
a vertex with incident set \(I\), every mirror through the vertex belongs to
the finite standard parabolic root system generated by the eight roots in
\(I\).  Exact reflection closure at all 41 vertices therefore exhausts every
touching mirror.

The sum of the 41 local mirror counts is the proved finite pre-deduplication
bound 1856.  Deduplication gives 320 restricted mirrors:

| contact-face dimension | mirrors |
|---:|---:|
| 0 | 88 |
| 1 | 43 |
| 2 | 26 |
| 3 | 17 |
| 4 | 96 |
| 5 | 25 |
| 6 | 14 |
| 7 (facets) | 11 |

Thus 309 mirrors merely touch a proper nonfacet face.  In signed-vector form,
the complete five-shell answer is

| field norm | signed vectors whose restricted hyperplane meets \(\overline P\) |
|---|---:|
| \(5-3\phi\) | 0 |
| \(4-2\phi\) | 640 |
| \(3-\phi\) | 0 |
| \(2\) | 640 |
| \(1+\phi\) | 0 |

The script materializes every coordinate vector.  The compact JSON archive
records the SHA-256 fingerprint
`815c0c4cad8e4ab9da32ef491c104ac71a781d4f2b4300f988ecb484506da69e`
of the 320 canonical norm-two representatives.

## Weyl-conjugate nef realization

Choose a point of \(P^\circ\), and let \(D\) be the chamber of the full K3
root arrangement in the positive cone of \(N\) which contains it.  Equality
of the restricted arrangements gives

\[
D\cap E=P.
\]

For a generic K3 surface with \(\operatorname{NS}(X)=N\), let \(A_X\) be its
ample chamber.  The K3 Weyl group is transitive on root chambers, so there is
\(w\in W(N)\) with \(wD=A_X\).  Conjugate the endomorphism structure,

\[
\iota_w(a)=w\iota(a)w^{-1},\qquad E_w=wE.
\]

Then

\[
A_X\cap E_w=wP,
\qquad
\operatorname{Nef}(X)\cap\overline{E_w^+}=w\overline P.
\]

This proves the existence of a marked K3 eigen-nef section isometric to the
compact Bugaenko Coxeter polytope.

The conjugation is compatible with the fixed primitive K3 embedding.  For a
K3 root \(d\in N\) and \(x\in N^\vee\),

\[
s_d(x)-x=(x,d)d\in N.
\]

Hence every root reflection, and therefore every \(w\in W(N)\), acts
trivially on the discriminant group \(A_N\).  The map
\(w\oplus\operatorname{id}_{N^\perp}\) preserves the unimodular gluing and
extends to an isometry of \(\Lambda_{K3}\) fixing the period space.  This is a
change of marking/NS endomorphism structure, not an assertion that \(w\) is
induced by an automorphism of \(X\).

Finally, the restricted reflection arrangement is invariant, but its chamber
is not.  If \(h(r,r)=2\), put \(u=\phi-1\).  The two trace roots \(r,ur\)
are orthogonal and their K3 root reflections multiply to the
\(\mathcal O_F\)-linear field reflection.  Explicitly, for
\(c=h(x,r)\),

\[
\operatorname{Tr}(\alpha c)
+u\operatorname{Tr}(\alpha u c)=-c;
\]

the certificate checks this on the basis \(c=1,\phi\).  Field reflections
preserve \(E\) and permute all restricted mirrors.  A facet reflection sends
\(P\) to the adjacent chamber; it does **not** preserve the nef section
setwise.  Fundamental-domain invariance of the arrangement, not setwise
invariance of one chamber, is the Coxeter statement.

## Reproduction

```sh
PYTHONPATH=src:. python3 scripts/bugaenko_k3_eigen_nef.py
PYTHONPATH=src:. pytest -q tests/test_bugaenko_k3_eigen_nef.py
```

The summary certificate is
`data/trace_lattices/bugaenko_k3_eigen_nef.json`.
