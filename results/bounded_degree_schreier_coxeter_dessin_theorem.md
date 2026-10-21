# Bounded-degree Schreier--Coxeter dessins and the dotted-height boundary

**Status (2026-08-09): proved arithmetic/group-theoretic reduction, with an
exact rational ideal-polygon counterfamily.** This closes the formerly vague
dotted-height step after one indispensable complexity parameter is added:
the degree of the chamber cover, equivalently relative covolume. It does not
bound that degree from Borcherds data, and therefore does not by itself give a
new dimension bound.

## 1. The finite dessin

Let \(P\subset\mathbb H^n\) be a compact or finite-volume arithmetic Coxeter
polytope and let \(W=W(P)\). Choose a maximal arithmetic reflection overgroup

\[
 W\leq W_0=W(Q),\qquad
 I=[W_0:W]=\frac{\operatorname{vol}P}{\operatorname{vol}Q}.
\]

Different maximal reflection overgroups can exist. The canonical numerical
complexity is

\[
 \deg_{\mathrm{Cox}}(P)=
 \min_{W\leq W_0\ {\rm maximal\ reflection}}[W_0:W].
\]

The minimum exists because it is over positive integers. Choose a minimizing
\(W_0\) when canonical degree matters; a tie between minimizing overgroups is
genuine noncanonicity and is retained in the dessin.

Such an overgroup exists: every discrete overgroup of the lattice \(W\) has
finite index, an ascending chain has bounded length by the hyperbolic-orbifold
lower-volume bound, and a maximal reflection overgroup remains arithmetic.
Agol--Belolipetsky--Storm--Whyte prove that only finitely many conjugacy
classes of maximal arithmetic reflection groups occur.

Write the finite Coxeter generating set of \(W_0\) as
\(S=\{s_1,\ldots,s_e\}\). The **Schreier--Coxeter dessin** of the based
inclusion \(W\leq W_0\) is

\[
 \mathfrak D(P/Q)=
 \bigl(Q,\Omega=W\backslash W_0,\omega=W,
       \sigma_i:Wx\longmapsto Wxs_i\bigr).
\]

It is a finite transitive permutation tuple. Every \(\sigma_i\) is an
involution and the Coxeter relations of \(W_0\) hold. This is an exact, not
metaphorical, analogue of the monodromy-permutation form of a dessin:

- the degree is \(I=|\Omega|\);
- a transposition of \(\sigma_i\) glues two adjacent \(Q\)-tiles;
- a fixed flag \(Wx=Wxs_i\) is a mirror flag because
  \(xs_ix^{-1}\in W\);
- the stabilizer of the marked sheet \(\omega\) is exactly \(W\).

The finite reflection gate is also exact. Reidemeister--Schreier gives a
finite word generating set for the stabilizer \(W\). For every fixed flag,
choose a word representative \(x\) and retain the reflection
\(xs_ix^{-1}\). A finite list of word expressions of the Schreier generators
in these reflections certifies

\[
 W=\left\langle
 W\cap\{xs_ix^{-1}:x\in W_0,\ s_i\in S\}
 \right\rangle .
\]

Conversely this equality says precisely that the stabilizer is a reflection
subgroup. Its geometric chamber, tiled by the \(I\) copies of \(Q\), recovers
\(P\). Thus a maximal base chamber, a finite transitive Coxeter permutation
tuple, and the finite reflection certificate determine \(P\) up to based
isometry. Dyer's canonical-generator theorem is an equivalent intrinsic
route to the Coxeter generating reflections.

For unbased classification, quotient by ambient isometries and wall
permutations. Different choices of a base tile inside the same \(P\) are
extra markings, not new unmarked polytopes.

This is an arithmetic reflection/coset dessin. It is **not** a claim that
\(W\backslash\mathcal D\) is a finite algebraic type-IV cover: \(W\) need
not be a lattice in the signature-\((2,n)\) group. The existing
determinant-line Shimura variety and marked Heegner data remain the algebraic
carrier.

## 2. Bounded degree gives an explicit all-embedding height bound

Fix an arithmetic realization of \(W_0\) over its totally real field \(k\),
including its coordinate fractional ideals and a common denominator \(D_Q\).
Equivalently, work in a fixed \(W_0\)-stable free \(\mathcal O_k\)-lattice
when one is available. Let \(a_i\) be the base roots and \(S_i\) the
reflection matrices. Put \(d=[k:\mathbb Q]\), choose a submultiplicative
matrix norm at every real embedding, and set

\[
 M_Q=\max\left(1,\|\sigma(S_i)\|:
        \sigma:k\hookrightarrow\mathbb R,\ 1\leq i\leq e\right).
\]

There is an explicit basis-dependent constant \(C_Q\geq1\), computed from
the base roots and bilinear matrix, such that

\[
 \left|\sigma\bigl((a_i,ga_j)\bigr)\right|
 \leq C_Q\|\sigma(g)\|
 \quad\text{for every }i,j,\sigma,g.
\]

> **Bounded-degree height theorem.** Every facet root of an index-\(I\)
> reflection-subgroup chamber \(P\supset Q\) is \(ga_i\), where the tile
> \(gQ\subset P\) can be chosen at gallery distance at most \(I-1\) from
> \(Q\). Hence every integral facet Gram entry satisfies, at every embedding,
>
> \[
> \left|\sigma\bigl((r_u,r_v)\bigr)\right|
> \leq B_Q(I):=C_QM_Q^{\,2I-2}.
> \]
>
> Facet-root coordinates are bounded by \(R_QM_Q^{I-1}\) for another
> explicit base constant \(R_Q\), and facet-reflection matrices are bounded
> by \(M_Q^{2I-1}\).

Indeed, the \(I\) base chambers inside the convex chamber \(P\) have a
connected gallery graph. A shortest spanning-tree path to any tile has
length at most \(I-1\). A facet of \(P\) contains a wall-piece of one such
tile, so its reflection is \(gS_ig^{-1}\) and its root is \(ga_i\).
For two facet roots, \(g^{-1}h\) has word length at most \(2I-2\), which
gives the displayed bound. The proof is identical for compact and
finite-volume noncompact chambers; compactness is invariant under this
finite cover.

The root norm and divisibility alphabet is automatically finite here: every
facet reflection of \(W\) is conjugate in \(W_0\) to a base simple
reflection. Prescribing a smaller alphabet is an additional finite filter.

After multiplying coordinates and normalized Gram entries by the retained
denominator \(D_Q\), choose the rational principal level
\(m\mathcal O_k\) with

\[
 m>2D_Q\max\{B_Q(I),R_QM_Q^{I-1}\}.
\]

The norm inequality for algebraic integers then makes reduction injective
on every allowed coordinate and Gram entry and separates distinct root
lines up to sign. Thus the **carrier-marked** finite dessin--the reduced
boundary roots together with their tiling/Schreier marking--recovers the
exact Gram. The unmarked image of \(W_0\) in a congruence quotient need not
recover a noncongruence subgroup \(W\). The literal Vinberg
compactness/finite-volume test remains terminal.
This supplies an explicit upper level bound for every bounded-degree class.
A smaller good ideal can work in an individual case.

## 3. Classification pushback

For a fixed maximal base \(Q\) and \(I\leq I_0\), only finitely many
transitive tuples of permutations of degree at most \(I_0\) satisfy its
Coxeter relations. The finite reflection gate selects exactly the reflection
subgroups. Therefore the corresponding unmarked \(P\)'s, modulo ambient
isometry and wall permutation, are finite and effectively enumerable from
\(Q\); their exact root systems, congruence levels,
determinant-line carriers, obstruction matrices, and Vinberg volume
certificates are computable.

An explicitly terminating implementation need not solve an abstract subgroup
equality problem: enumerate the finite connected subsets of at most \(I_0\)
chambers in the gallery ball of radius \(I_0-1\) about \(Q\), then test exact
convexity, Coxeter boundary angles, and the complete Vinberg volume criterion.
Every chamber in the theorem occurs in this list.

Globally, impose \(\operatorname{vol}P\leq V\). The ABSW finiteness theorem
gives finitely many maximal arithmetic bases \(Q\), and for each of them

\[
 I\leq\left\lfloor
 \frac{V}{\operatorname{vol}Q}
 \right\rfloor .
\]

Consequently there are only finitely many arithmetic Coxeter \(P\)'s of
bounded covolume, and every one has a finite exact Schreier--Coxeter dessin
and an existential uniform congruence-height bound. This global statement
is not yet an effective census: ABSW proves the maximal-base list is finite
but does not supply the compiled list and constants used here. With an
explicit base list it becomes a terminating classification.

The field \(k\), its discriminant, dimension \(n\), compact/noncompact type,
and root alphabet are labels of the base and useful filters. They do **not**
replace \(I\) or \(V\). The field is unchanged under finite index. Known
field-degree or discriminant estimates for maximal or congruence reflection
groups can shrink the finite base list, but do not control the subgroup
index.

## 4. Exact no-go: \(n,k\), the lattice and root alphabet do not bound height

Use the hyperboloid with form

\[
 \langle x,y\rangle=-x_0y_0+x_1y_1+x_2y_2
\]

and the arithmetic ideal triangle
\(\Delta_0=(\infty,0,1)\). Its reflection group lies in
\(O(I_{2,1})\). For \(N\geq1\), set

\[
 P_N=\bigcup_{j=0}^{N-1}(\infty,j,j+1).
\]

This is the convex ideal \((N+2)\)-gon bounded by \(x=0\), \(x=N\), and the
semicircles with endpoints \(j,j+1\). It is a union of \(N\) base chambers,
so

\[
 [W(\Delta_0):W(P_N)]=N,\qquad
 \operatorname{area}P_N=N\pi.
\]

The semicircle wall \([j,j+1]\) has primitive norm-one integral root

\[
 n_j=\bigl(-1-j(j+1),-(2j+1),1-j(j+1)\bigr).
\]

Under the standard embedding

\[
 X(x+iy)=\frac{1}{2y}
 (x^2+y^2+1,2x,x^2+y^2-1),
\]

the equation \(\langle X,n_j\rangle=0\) is
\(x^2+y^2-(2j+1)x+j(j+1)=0\). Since a norm-one integral root acts by
\(v\mapsto v-2\langle v,r\rangle r\), all these wall reflections preserve
\(I_{2,1}\). The vertical walls have primitive norm-one roots
\((0,1,0)\) and \(-(N,1,N)\). Direct calculation gives

\[
 \langle n_j,n_j\rangle=1,\qquad
 \langle n_0,n_j\rangle=1-2j^2.
\]

For \(j\geq2\), these are ultraparallel walls and
\(\cosh d([0,1],[j,j+1])=2j^2-1\). Hence dimension \(2\), field
\(\mathbb Q\), the odd unimodular lattice \(I_{2,1}\), finite-volume
noncompact type, and the primitive norm/divisibility alphabet
\(\{(1,1)\}\) all stay fixed while facet number, volume and dotted Gram
height are unbounded.

There is also no universal fixed congruence level. For every \(m\geq2\),

\[
 n_{j+m}\equiv n_j\pmod m.
\]

Thus level \(m\) fails to separate facet sheets as soon as both occur. For
an arbitrary fixed ideal \(\mathfrak a\subset\mathcal O_k\), take a positive
integer \(m\in\mathfrak a\cap\mathbb Z\) and the same periodicity applies.

This is the sharp no-go boundary: no height, facet, volume,
number-of-marked-walls, or single-level faithful-dessin bound can depend
only on \(n,k,\operatorname{disc}k\), compactness type and a root alphabet.
A surface-degree claim would require specifying which marked curves and
monodromy complexity its Hilbert polynomial controls. Allcock's redoubling
theorem supplies the same qualitative obstruction in fixed arithmetic
commensurability classes through much higher dimensions, including compact
sequences; the elementary \(P_N\) family is the exact-coordinate certificate.

## 5. Consequence for Borcherds dessins

Unmarked ambient lattice/Borcherds data cannot distinguish arbitrary
finite-index reflection subgroups of the same arithmetic group. The missing
finite datum is the chamber-cover monodromy \(\mathfrak D(P/Q)\). Once its
degree \(I\) is bounded, the theorem above bounds dotted height and a faithful
congruence level, so the determinant-line construction turns each of
finitely many candidates into a finite obstruction problem.

What remains open is the genuinely automorphic implication

\[
 \boxed{\text{bounded \(P\)-specific Borcherds complexity}
 \Longrightarrow I\leq F(n,k,L,\text{complexity}).}
\]

No such inequality is currently proved. It cannot use only the unmarked
ambient lattice, because the family above keeps that lattice fixed. A viable
statement must use a quantity that grows under chamber doubling: cover
degree, marked-sheet degree, automorphic degree on a sheet-separating cover,
or equivalent monodromy complexity. Proving this inequality would complete
the route:

    bounded P-specific Borcherds degree
      -> bounded Schreier degree I
      -> explicit dotted height and congruence level
      -> finite permutation-dessin census
      -> exact Gram lift and Vinberg test
      -> examples or a nonexistence/facet theorem

## 6. Frozen exact certificate and sources

The ideal-polygon identities and congruence collision are replayed by
[the exact certificate script](../scripts/certify_ideal_triangle_coxeter_dessin_no_go.py)
and its
[frozen JSON archive](../data/enumerations/ideal_triangle_coxeter_dessin_no_go.json).

Primary sources:

- Agol--Belolipetsky--Storm--Whyte,
  [*Finiteness of arithmetic hyperbolic reflection groups*](https://arxiv.org/abs/math/0612132).
- Allcock,
  [*Infinitely many hyperbolic Coxeter groups through dimension 19*](https://arxiv.org/abs/0903.0138).
- Felikson--Tumarkin,
  [*On finite index reflection subgroups of discrete reflection groups*](https://arxiv.org/abs/math/0305093).
- Dyer, *Reflection subgroups of Coxeter systems*, Journal of Algebra 135
  (1990), 57--73, DOI 10.1016/0021-8693(90)90149-I.
