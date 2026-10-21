# Compact right-angled K3 hexagons have a finite field set

**Status (2026-08-11): proved for every full-trace complex-K3 field
degree.**  Let a rank-three module over a totally real field return an exact
compact right-angled Coxeter hexagon in \(\mathbf H^2\).  The full trace has
rank \(3d\), so

\[
                              d=[k:\mathbf Q]\le6.
\tag{0.1}
\]

Rank-three Gram closure is the right-hexagon cosine recurrence.  It implies
that among the six second-neighbour dotted pairs there is one whose shifted
squared ratio satisfies

\[
                                  0<x_i<16.
\tag{0.2}
\]

Its pair cyclic Cartan product

\[
                                  p_i=4(x_i+1)
\tag{0.3}
\]

is an algebraic integer with distinguished value in \((4,68)\), while all
its other conjugates lie in \((0,4)\).  This unique-large-conjugate pattern
forces \(p_i\) itself to generate \(k\).  Consequently,

\[
                    |\operatorname {Disc}(k)|\le68^{30}.
\tag{0.4}
\]

Hermite--Minkowski gives a finite field set.  No imposed height, power
coupling, Coxeter classification, reflective-lattice classification,
saturation, or maximality enters the bound.

## 1. Exact carrier package

Let \(k\) be totally real of degree \(d\), let \(\mathcal O\subset k\) be a
returned order, and let \(M\) be a finite projective rank-three
\(\mathcal O\)-module with an admissible symmetric \(k\)-form \(h\).  At one
distinguished embedding \(\sigma_0\), the form has signature \((2,1)\); at
every other embedding it is positive definite.  Assume:

1. a minimal complex projective K3 surface \(X\) has Picard lattice equal
   to the primitive full trace of \((M,h)\), and the returned order acts
   self-adjointly by divisor-product correspondences;
2. the prescribed projectivized eigen-nef chamber is exactly a compact
   right-angled Coxeter hexagon in \(\mathbf H^2\);
3. all six facet reflections preserve \(M\); and
4. \(k\) is the minimal cyclic Cartan ground field of the returned
   reflection packet.

Complete saturation and declared-host reflection maximality may be retained
but are not used.  The primitive full trace gives
\(3d=\rho(X)\le20\), proving (0.1).  If \(d=1\), the field is already
\(\mathbf Q\); below one may assume \(2\le d\le6\).

## 2. Six second-neighbour coordinates

Label the sides cyclically by \(L_0,\ldots,L_5\).  Consecutive sides meet
orthogonally.  The sides \(L_i,L_{i+2}\) are ultraparallel; let

\[
 A_i={h(e_i,e_{i+2})^2\over h(e_i,e_i)h(e_{i+2},e_{i+2})},
 \qquad x_i=A_i-1>0.
\tag{2.1}
\]

Geometrically, \(A_i=\cosh^2\ell_i\) and
\(x_i=\sinh^2\ell_i\), where \(\ell_i\) is the intervening side length.
Only the squared ratios \(A_i\) are asserted to lie in \(k\); the positive
square roots below are distinguished-place real numbers.

Choose the three alternating variables \(x_0,x_2,x_4\), and put

\[
 A=x_0+1,
 \qquad C=x_2+1,
 \qquad E=x_4+1,
 \qquad T=\sqrt{ACE}.
\tag{2.2}
\]

The rank-three Gram calculation for the alternating normals gives the
right-hexagon closure

\[
\begin{aligned}
 S&=A+C+E-1+2T,\\
 x_1&={S\over x_0x_2},
 &x_3&={S\over x_2x_4},
 &x_5&={S\over x_4x_0}.
\end{aligned}
\tag{2.3}
\]

For completeness, if the alternating normalized correlation magnitudes are
\(u=\sqrt A\), \(v=\sqrt C\), \(w=\sqrt E\), their Gram matrix is

\[
 \begin{pmatrix}
 1&-u&-w\\
 -u&1&-v\\
 -w&-v&1
 \end{pmatrix}.
\tag{2.4}
\]

The normal to the first adjacent pair and the normal to the next adjacent
pair are dual cofactor vectors.  Their squared normalized correlation is

\[
 { (uv+w)^2\over(u^2-1)(v^2-1)}
 ={AC+E+2T\over x_0x_2}=x_1+1,
\tag{2.5}
\]

which proves the first identity in (2.3); rotation proves the other two.
This is a rank-three cofactor calculation, not an imported hexagon
classification.

Four-consecutive Gram minors give the compatible relations

\[
                         x_i x_{i+1}=x_{i+3}x_{i+4}
                         \qquad(i\bmod6).
\tag{2.6}
\]

## 3. A uniformly short dotted pair

If one of \(x_0,x_2,x_4\) is below \(16\), equation (0.2) already holds.
Otherwise order the alternating variables as

\[
                     x\le y\le z,
                     \qquad x,y,z\ge16.
\]

The odd variable whose denominator in (2.3) is \(yz\) satisfies, with
\(C=y+1\) and \(E=z+1\),

\[
\begin{aligned}
 x_{\mathrm{odd}}
 &={A+C+E-1+2\sqrt{ACE}\over yz}\\
 &<{2C+E+2C\sqrt E\over yz}\\
 &\le {51+34\sqrt{17}\over256}\\
 &<1.
\end{aligned}
\tag{3.1}
\]

Here \(A\le C\le E\), while the functions
\((t+1)/t\) and \(\sqrt{t+1}/t\) decrease for \(t>0\).  The last inequality
uses \(\sqrt{17}<5\), so its numerator is below \(221<256\).
Thus in every case some one of the six \(x_i\) is below \(16\), proving
(0.2).

## 4. The dotted pair is a primitive field generator

Fix the index supplied by Section 3.  The pair cycle

\[
 p_i=4A_i=4(x_i+1)
\tag{4.1}
\]

belongs to \(\mathcal O\): it differs by an integer from the module trace
of the product of the two integral facet reflections.  At the distinguished
place, compact ultraparallelity and (0.2) give

\[
                              4<(p_i)_0<68.
\tag{4.2}
\]

At every nondistinguished embedding \(\sigma\), positive-definite
Cauchy--Schwarz gives

\[
                              0<\sigma(p_i)<4.
\tag{4.3}
\]

The lower endpoint is strict because the original dotted pairing is a
nonzero field element and embeddings are injective.

### Lemma 4.1 (unique-large-conjugate primitivity)

Let \(k\) be totally real, and let \(p\in k\) have one real embedding with
value greater than \(4\) and every other embedding with value below \(4\).
Then \(\mathbf Q(p)=k\).

Indeed, put \(F=\mathbf Q(p)\).  If \([k:F]>1\), the restriction of the
distinguished embedding to \(F\) has \([k:F]\) distinct extensions to
\(k\).  All are real because \(k\) is totally real, and all take the same
value \((p_i)_0>4\) on \(p\).  At least one of them is nondistinguished,
contradicting (4.3).  Hence \(F=k\).

This argument also shows why a gratuitous proper totally real field
extension cannot hide above a dotted pair in a one-Lorentz-place admissible
carrier.

## 5. Discriminant and theorem

The algebraic integer \(p_i\) is primitive by Lemma 4.1, and all its
conjugates lie in \((0,68)\).  Thus every conjugate difference has magnitude
below \(68\), and

\[
\begin{aligned}
 |\operatorname {Disc}(k)|
 &\le|\operatorname {Disc}(\mathbf Z[p_i])|\\
 &<68^{d(d-1)}\\
 &\le68^{30}.
\end{aligned}
\tag{5.1}
\]

Hermite--Minkowski proves (0.4) and the finite field set.

### Theorem 5.1 (compact right-hexagon field theorem)

Every exact full-trace complex-K3 carrier satisfying Section 1 and having
a compact right-angled hexagon as its prescribed eigen-nef chamber belongs
to the finite set of totally real fields bounded by (0.4).

## 6. Scope and the higher-facet frontier

The theorem does not cover non-right hexagons, noncompact hexagons, or
arbitrary chambers with six or more facets.  It neither constructs the K3
carrier nor proves essential surjectivity.  Saturation and maximality are
not needed.

Lemma 4.1 gives a useful general reduction: for any exact admissible
one-Lorentz-place carrier, a dotted pair cycle is already a primitive field
generator.  Therefore an arbitrary compact polygon would have a finite
field set once OUR geometry supplies a uniform upper bound for at least one
dotted pair.  Proving that short-dotted-pair statement beyond the present
right-hexagon incidence is the remaining geometric obligation.
