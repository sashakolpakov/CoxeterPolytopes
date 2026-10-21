# Rank-four K3 correspondence packages carry infinitely many real quadratic fields

**Status (2026-08-10): proved scoped no-go.**  A fixed Neron--Severi rank,
minimality, a primitive full trace lattice, an integral self-adjoint field
action by algebraic correspondence, and exact equality with a complete
restricted K3 eigen-nef component do **not** force a finite set of totally
real fields.

There is an explicit family indexed by every integer \(a\ge2\):

\[
 k_a=\mathbf Q\bigl(\sqrt{a^2-1}\bigr),
 \qquad \rho(X_a)=4,                                   \tag{0.1}
\]

and infinitely many of the fields \(k_a\) are distinct.  Each package is a
minimal projective K3 surface with exact Picard lattice, a primitive full
rank-four restriction-of-scalars block, and a divisor-product
correspondence whose action on that Picard block generates \(k_a\).

This theorem deliberately stops short of the strongest field-set gate.  It
does **not** prove that the restricted eigen-nef component is finite-volume,
Coxeter, reflection-saturated, reflection-maximal, or equal to a prescribed
input chamber.  Thus literal reflection-maximal same-\(P\) eigen-nef return
could still impose additional arithmetic bounds.  What is ruled out is any
argument using bounded rank, minimal K3 geometry, primitive trace embedding,
or correspondence existence alone.

## 1. The integral family

Fix \(a\ge2\) and put

\[
 p_a(x)=x^2-2ax+1,
 \qquad
 t_a=a+\sqrt{a^2-1},
 \qquad t_a'=a-\sqrt{a^2-1}.                           \tag{1.1}
\]

Then \(t_a>1\), \(0<t_a'<1\), and \(p_a\) is irreducible: the integer
\(a^2-1\) lies strictly between \((a-1)^2\) and \(a^2\).  Let

\[
 C_a=\begin{pmatrix}0&-1\\1&2a\end{pmatrix},
 \quad
 B_0=\begin{pmatrix}-1&0\\0&1\end{pmatrix},
 \quad
 B_a=\begin{pmatrix}1-2a&-1\\-1&-1\end{pmatrix}.     \tag{1.2}
\]

Define the even integral rank-four lattice

\[
 L_a=\bigl(\mathbf Z^4,G_a\bigr),
 \qquad
 G_a=2B_0\oplus2B_a,                                  \tag{1.3}
\]

and the integral operator

\[
                         \tau_a=C_a\oplus C_a.         \tag{1.4}
\]

### Proposition 1.1 (full trace and signatures)

The package \((L_a,G_a,\tau_a)\) has the following exact properties.

1. \(G_a\) is even of signature \((1,3)\), and
   \[
                              \det G_a=-32(a-1).        \tag{1.5}
   \]
2. \(\tau_a\) is self-adjoint and satisfies
   \[
                         \tau_a^2-2a\tau_a+I=0.        \tag{1.6}
   \]
   Hence \(\mathbf Q[\tau_a]\simeq k_a\).
3. With \(A_a=\mathbf Z[t_a]\), the lattice \(L_a\) is a free
   \(A_a\)-module of rank two.  Thus its underlying rational vector space is
   the full restriction of scalars of \(k_a^2\), not one real eigenspace.
   Moreover the self-adjoint rational form is the trace of a unique
   \(k_a\)-bilinear symmetric form, so this is a literal full trace package.
4. The \(t_a\)-eigenspace has signature \((1,1)\), while the
   \(t_a'\)-eigenspace is negative definite.

### Proof

The block \(B_0\) has signature \((1,1)\).  For \(B_a\),

\[
                  1-2a<0,
          \qquad \det B_a=2(a-1)>0,                   \tag{1.7}
\]

so it is negative definite.  This proves the signature.  Scaling the two
binary blocks gives

\[
       \det(2B_0)=-4,
       \qquad \det(2B_a)=8(a-1),                      \tag{1.8}
\]

which proves (1.5).

Direct multiplication gives

\[
 C_a^tB_0=B_0C_a=
       \begin{pmatrix}0&1\\1&2a\end{pmatrix},
 \qquad
 C_a^tB_a=B_aC_a=
       \begin{pmatrix}-1&-1\\-1&1-2a\end{pmatrix}.  \tag{1.9}
\]

Thus \(\tau_a\) is self-adjoint.  Its two blocks are the regular action of
\(t_a\) on \(A_a\) in the basis \((1,t_a)\), and (1.6) follows.  The two
blocks make \(L_a\simeq A_a^2\).  More intrinsically, for \(x,y\in
L_a\otimes\mathbf Q\), the functional
\(c\mapsto (cx).y\) on \(k_a\) has the form
\(\operatorname {Tr}_{k_a/\mathbf Q}(c\,h(x,y))\) for a unique \(h(x,y)\in
k_a\).  Self-adjointness of the generator makes \(h\) \(k_a\)-bilinear.
This supplies the claimed trace form without discarding the conjugate real
place.

For an eigenvalue \(\lambda\in\{t_a,t_a'\}\), the vector
\((1,-\lambda)\) is a \(C_a\)-eigenvector.  Its \(B_0\)-square is

\[
                              \lambda^2-1.              \tag{1.10}
\]

This is positive for \(t_a\) and negative for \(t_a'\).  The corresponding
line in the \(B_a\)-block is always negative because \(B_a\) is negative
definite.  The claimed two real-place signatures follow. \(\square\)

## 2. Explicit primitive embedding in the K3 lattice

Let

\[
 \Lambda_{K3}=U_1\oplus U_2\oplus U_3
                 \oplus E_8(-1)\oplus E_8(-1),        \tag{2.1}
\]

where \(e_i,f_i\) is the standard basis of \(U_i\), so
\(e_i^2=f_i^2=0\) and \(e_i.f_i=1\).  Choose a primitive root
\(v\in E_8(-1)\) with \(v^2=-2\).

In the block basis \((x_1,x_2,y_1,y_2)\) of (1.3), define

\[
\begin{aligned}
 x_1&\longmapsto e_1-f_1,\\
 x_2&\longmapsto e_2+f_2,\\
 y_1&\longmapsto e_3+(2-2a)f_3+v,\\
 y_2&\longmapsto v.
\end{aligned}                                          \tag{2.2}
\]

### Proposition 2.1

Equation (2.2) is a primitive isometric embedding

\[
                            L_a\hookrightarrow\Lambda_{K3}.     \tag{2.3}
\]

### Proof

The four image squares are

\[
                       -2,\quad2,\quad2-4a,\quad-2,    \tag{2.4}
\]

and the only nonzero off-diagonal pairing is

\[
                    \bigl(e_3+(2-2a)f_3+v\bigr).v=-2. \tag{2.5}
\]

This is exactly \(2B_0\oplus2B_a\).

Primitivity is visible without a discriminant-form theorem.  Replace the
third image vector by its difference with the fourth.  The resulting four
vectors have coefficient-one pivots on

\[
                         e_1,\quad e_2,\quad e_3,\quad v.         \tag{2.6}
\]

Because the root \(v\) is primitive and extends to an integral basis of
\(E_8(-1)\), the corresponding \(4\times4\) minor has determinant one.
Thus the quotient of \(\Lambda_{K3}\) by the image is torsion-free. \(\square\)

No lattice-embedding classification is used here.

## 3. Minimal projective K3 carriers and correspondences

The orthogonal complement of (2.3) has signature \((2,16)\).  Choose a
period in its type-IV domain outside the countable union of rational
Noether--Lefschetz hyperplanes.  K3 period surjectivity gives a marked K3
surface \(X_a\) with

\[
                         \operatorname {NS}(X_a)\simeq L_a.      \tag{3.1}
\]

The class \(x_2\) has square two, so \(X_a\) is projective.  It is minimal
and \(K_{X_a}=0\).

The integral self-adjoint operator \(\tau_a\) is induced by a
rational algebraic correspondence supported on products of divisors.  (The operator
is integral on \(\operatorname {NS}\); the cycle need only have rational
coefficients.)  Indeed, choose
a basis \(D_i\) of \(\operatorname {NS}(X_a)_{\mathbf Q}\) and its
intersection-dual basis \(D_i^\vee\).  Then

\[
              Z_{\tau_a}=\sum_i D_i^\vee\times\tau_a(D_i)
                  \in CH^2(X_a\times X_a)_{\mathbf Q}            \tag{3.2}
\]

acts on \(\operatorname {NS}(X_a)_{\mathbf Q}\) as \(\tau_a\), because

\[
 (D\times E)_*(\xi)=(\xi.D)E.                         \tag{3.3}
\]

Thus \(\mathbf Q[Z_{\tau_a}]\) supplies the complete field action on the
full trace block.

### Exact restricted eigen-nef component

Let \(E_a\subset L_a\otimes\mathbf R\) be the \(t_a\)-eigenspace.  It has
signature \((1,1)\).  No nonzero rational vector lies in the conjugate
eigenspace: an irrational eigenvalue cannot have a nonzero rational
eigenvector.  In particular no integral K3 root is perpendicular to all of
\(E_a\), so every restricted root mirror is a proper hyperplane.  Choose a
positive regular point \(p\in E_a\) avoiding their countable union, and let
\(C_a^{\rm res}\) be its component in the complete restricted root
arrangement.  Let \(D\) be the ambient K3 root chamber containing \(p\).
The root arrangement is locally finite in the positive cone, and equality
of sign patterns gives

\[
                              D\cap E_a=C_a^{\rm res}.             \tag{3.4}
\]

The K3 Weyl group is transitive on ambient root chambers.  Choose a Weyl
element \(w\) carrying \(D\) to the ample chamber and conjugate the field
action by \(w\).  With \(E_{a,w}=wE_a\), convexity and the regular point in
the relative interior give

\[
 \boxed{
 \mathbf P\bigl(\operatorname {Nef}(X_a)\cap E_{a,w}^+\bigr)
          =\mathbf P\bigl(w\overline{C_a^{\rm res}}\bigr).}      \tag{3.5}
\]

This is exact equality with a component of the **complete** restricted K3
root arrangement.  It is not the unproved stronger assertion that this
component returns a prescribed arithmetic Coxeter chamber or is maximal.

## 4. Infinitely many fields

### Proposition 4.1

The set \(\{k_a:a\ge2\}\) contains infinitely many distinct real quadratic
fields.

### Proof

Use a finite-union contradiction.  Suppose only the squarefree fields
\(\mathbf Q(\sqrt{d_1}),\ldots,
\mathbf Q(\sqrt{d_s})\) occurred.  For every \(a\), the equality
\(k_a=\mathbf Q(\sqrt{d_j})\) would give an integer \(b\ge1\) with

\[
                              a^2-d_jb^2=1.             \tag{4.1}
\]

For a fixed \(d_j\), its positive Pell solutions are powers of a fundamental
unit.  Their first coordinates grow exponentially, so only
\(O_{d_j}(\log X)\) solutions have \(a\le X\).  A finite union of the fields
would therefore account for only \(O(\log X)\) integers \(a\le X\), whereas
the family contains all \(a=2,\ldots,X\).  This contradiction proves the
claim. \(\square\)

Combining Propositions 1.1, 2.1, and 4.1 proves the announced no-go.

## 5. What extra invariant makes the field image finite?

The rank inequality gives \([k:\mathbf Q]\le4/2=2\) in this family, yet the
field set is infinite.  Three exact repairs are available.

1. **Field discriminant.**  A uniform bound on \(|\operatorname {Disc}k|\),
   together with bounded degree, gives finitely many fields by
   Hermite--Minkowski.  This is the weakest direct scalar repair.
2. **Primitive spectral envelope.**  If one integral primitive generator
   \(\theta\) of \(k\) has degree at most \(R\) and every conjugate satisfies
   \(|\sigma(\theta)|\le H\), then its minimal-polynomial coefficients obey
   \[
          |c_{r-j}|\le {r\choose j}H^j.                \tag{5.1}
   \]
   Only finitely many polynomials and fields remain.
3. **Bounded integral operator packet.**  A uniform entry bound for a field-
   generating correspondence matrix in a bounded integral polarized lattice
   basis leaves only finitely many matrices and characteristic fields.  An
   operator norm without a bounded lattice/polarization normalization is
   basis-dependent and is not enough.

The present family escapes all three repairs:

\[
 \operatorname {disc}(p_a)=4(a^2-1),
 \qquad \det L_a=-32(a-1),
 \qquad \|C_a\|_{\max}=2a.                            \tag{5.2}
\]

A bound on the trace-lattice discriminant alone implies a field-discriminant
bound only after controlling the order index, codifferent trace scalar, and
determinant factor in the restriction-of-scalars formula.

## 6. Hostile scope boundary

The theorem proves:

```text
fixed NS rank 4
+ primitive full trace block
+ minimal projective K3 carrier
+ integral self-adjoint NS operator induced by a rational divisor correspondence
+ exact complete restricted-root eigen-nef component
does not imply a finite totally-real field set.
```

It does **not** prove any of the following:

- that the restricted component is a Coxeter polytope;
- that it is finite-sided, compact, or finite-volume;
- that its facet group is saturated in \(O(L_a)\);
- that it is reflection-maximal;
- that it returns a prescribed input chamber; or
- that exact reflection-maximal same-\(P\) carriers admit infinitely many
  fields.

Consequently the strongest exact-maximal field-set gate remains open.  Any
finite-field theorem there must prove that maximal same-\(P\) return supplies
one of the additional integral bounds in Section 5; bounded Picard rank and
correspondence existence do not supply it automatically.

The exact finite replay is supplied by the
[certifier](../scripts/certify_rank4_k3_correspondence_infinite_field_no_go.py),
[archive](../data/enumerations/rank4_k3_correspondence_infinite_field_no_go.json),
[tests](../tests/test_rank4_k3_correspondence_infinite_field_no_go.py), and
[hostile scope audit](rank4_k3_correspondence_infinite_field_no_go_hostile_scope_audit.md).
