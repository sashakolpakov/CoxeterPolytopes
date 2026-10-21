# Quarter-integral K3 right pentagons have a finite field set

**Status (2026-08-11): proved general compact-right-pentagon field
theorem.**  The power relation in the earlier
[scoped obstruction](k3_n2_power_coupled_right_pentagon_field_obstruction.md)
is unnecessary.  For an exact rank-six full-trace K3 return whose
prescribed eigen-nef chamber is a compact right-angled pentagon, the five
integral facet reflections alone make every shifted squared wall ratio
quarter-integral.  The five rank-three Gram identities form a cyclic
\(A_2\) recurrence.  A norm lower bound at the definite place then gives

\[
             {1\over16}<x_i<288\qquad(0\le i<5).
\tag{0.1}
\]

In particular, the moving closure denominator satisfies

\[
                         0<AC-A-C<288.
\tag{0.2}
\]

Every pair cyclic product already generates the returned real quadratic
field, and

\[
                   |\operatorname {Disc}(k)|\le1,343,277.
\tag{0.3}
\]

This is a finite field set with no imposed height and no relation between
\(A\) and \(C\).  Exact same-chamber return, saturation, and relative
reflection maximality are retained, although pair-reflection integrality
and right-pentagon closure prove the bound before maximality is used.

If the two initial squared ratios \(A,C\) themselves are algebraic
integers, the sharper bound

\[
                       |\operatorname {Disc}(k)|\le1152
\tag{0.4}
\]

holds, leaving 336 candidate quadratic squareclasses.  These are candidates,
not asserted realized K3 fields.

## 1. Exact geometric carrier package

Let \(k\) be a real quadratic field, \(\mathcal O\subset k\) a returned
order, and \(M\) a finite projective rank-three \(\mathcal O\)-module with
an admissible symmetric \(k\)-form \(h\).  At the distinguished embedding
\(\sigma_0\), the form has signature \((2,1)\); at the other embedding
\(\sigma_1\), it is positive definite.  Assume all of the following.

1. A minimal complex projective K3 surface \(X\) has rank-six Picard
   lattice equal to the primitive full trace of \((M,h)\).  The returned
   order acts self-adjointly through divisor-product correspondences.
2. The projectivized prescribed eigenslice of \(\operatorname {Nef}(X)\)
   is exactly a compact right-angled pentagon
   \(P\subset\mathbf H^2\), not merely a five-wall subsystem.
3. Its five facet reflections preserve \(M\).  The restricted K3 mirror
   packet is complete and saturated, and its group is reflection-maximal
   among integral \(k\)-linear reflections of the declared order host.

The projective nef dimension is \(\operatorname {rank}_k(M)-1=2\).  No
classification of reflective lattices is imported.

## 2. Five normalized pair coordinates

Choose cyclically ordered facet normals
\(e_0,\ldots,e_4\in M\otimes_{\mathcal O}k\).  Adjacent normals are
orthogonal.  For the five nonadjacent pairs, put

\[
\begin{array}{lll}
 A={h(e_0,e_2)^2\over h(e_0,e_0)h(e_2,e_2)},
 & C={h(e_1,e_3)^2\over h(e_1,e_1)h(e_3,e_3)},
 & B={h(e_0,e_3)^2\over h(e_0,e_0)h(e_3,e_3)},\\[3mm]
 D={h(e_1,e_4)^2\over h(e_1,e_1)h(e_4,e_4)},
 & E={h(e_2,e_4)^2\over h(e_2,e_2)h(e_4,e_4)}.
\end{array}
\tag{2.1}
\]

These invariant squared ratios lie in \(k\); no square root of a facet
norm is required to lie in \(k\).  At \(\sigma_0\), nonadjacent sides of a
compact right pentagon are ultraparallel, so

\[
                         A_0,B_0,C_0,D_0,E_0>1.
\tag{2.2}
\]

At \(\sigma_1\), strict positive-definite Cauchy--Schwarz gives

\[
                    0<A_1,B_1,C_1,D_1,E_1<1.
\tag{2.3}
\]

The lower endpoints are excluded because each nonadjacent pairing is a
nonzero element of \(k\) at \(\sigma_0\), hence remains nonzero under the
injective embedding \(\sigma_1\).  Equality at the upper endpoint would
make two normals dependent at \(\sigma_1\); their nonzero wedge over \(k\)
again excludes this.

## 3. All five Gram identities are the \(A_2\) pentagon

Every four facet normals lie in a three-dimensional space.  The five
principal four-by-four Gram determinants therefore vanish.  In the
notation (2.1), their equations are

\[
\begin{aligned}
 CE-C-D-E+1&=0,\\
 -A+BE-B-E+1&=0,\\
 BD-B-C-D+1&=0,\\
 AD-A-D-E+1&=0,\\
 AC-A-B-C+1&=0.
\end{aligned}
\tag{3.1}
\]

Define the shifted coordinates

\[
 x_0=A-1,\quad x_1=B-1,\quad x_2=C-1,
 \quad x_3=D-1,\quad x_4=E-1.
\tag{3.2}
\]

Equations (3.1) become the cyclic recurrence

\[
                  x_{i-1}x_{i+1}=x_i+1
                  \qquad(i\bmod5).
\tag{3.3}
\]

In particular, writing \(r=x_0=A-1\) and \(s=x_2=C-1\),

\[
\begin{aligned}
 x_1&=rs-1=AC-A-C,\\
 x_3&={s+1\over rs-1}=D-1,\\
 x_4&={r+1\over rs-1}=E-1.
\end{aligned}
\tag{3.4}
\]

Thus (3.3) uses all five closure identities, not only the two rational
formulas in (3.4).

## 4. Reflection integrality gives the quarter lattice

For a nonadjacent facet pair \(i,j\), its cyclic Cartan product is

\[
 a_{ij}a_{ji}
 ={4h(e_i,e_j)^2\over h(e_i,e_i)h(e_j,e_j)}.
\tag{4.1}
\]

It lies in \(\mathcal O\).  Indeed,
\(\operatorname {tr}_k(s_is_j)\) is the module trace of an endomorphism of
the finite projective \(\mathcal O\)-module and differs from (4.1) by an
integer.  Applying this to all five nonadjacent pairs gives

\[
                         4(x_i+1)\in\mathcal O,
             \qquad y_i:=4x_i\in\mathcal O.
\tag{4.2}
\]

No claim that \(x_i\in\mathcal O\) is made.  Equation (4.2), rather than
unjustified division by four, is the exact automatic integrality statement.

By (2.2)--(2.3),

\[
                 (y_i)_0>0,\qquad -4<(y_i)_1<0.
\tag{4.3}
\]

The algebraic integer \(y_i\) has negative nonzero rational norm.  Hence

\[
 1\le|N_{k/\mathbf Q}(y_i)|
   =(y_i)_0|(y_i)_1|<4(y_i)_0,
\]

which proves the uniform lower bound

\[
                              (x_i)_0>{1\over16}.
\tag{4.4}
\]

## 5. The cyclic upper bound

Fix an index \(i\) and put \(z=x_{i+2}\).  Two consecutive instances of
(3.3) give

\[
 x_{i+1}=x_i z-1,
 \qquad
 x_{i+3}={z+1\over x_i z-1}.
\tag{5.1}
\]

All distinguished coordinates are positive.  Combining
\((x_{i+3})_0>1/16\) with (5.1) yields

\[
 {1\over16}
 <{z+1\over x_i z-1},
 \qquad
 x_i<16+{17\over z}.
\tag{5.2}
\]

Since \(z>1/16\),

\[
                                  x_i<288.
\tag{5.3}
\]

The index was arbitrary, so (4.4) and (5.3) prove (0.1) for every one of
the five cyclic coordinates.  In particular
\(x_1=AC-A-C\) proves (0.2).  This is the missing moving-denominator bound.

The same estimate bounds every pair cyclic product at the distinguished
place:

\[
                    4<4(x_i+1)<1156.
\tag{5.4}
\]

The directed length-five Cartan cycle is bounded as well, since its square
is

\[
                    c_5^2=4^5ABCDE,
\]

and each of the five factors is below \(289\).  No separate mixed-height
premise enters.

## 6. The finite field theorem

Take any pair coordinate and write

\[
                           \alpha=4A=4(x_0+1)\in\mathcal O.
\tag{6.1}
\]

Equations (2.2)--(2.3) and (5.3) give

\[
                         4<\alpha_0<1156,
                  \qquad 0<\alpha_1<4.
\tag{6.2}
\]

In particular \(\alpha\notin\mathbf Q\): a rational element would have
equal values at the two embeddings.  Since \(k\) is quadratic,
\(k=\mathbf Q(\alpha)\).  Its integral trace and norm satisfy

\[
                 5\le T=\alpha_0+\alpha_1\le1159,
        \qquad 1\le N=\alpha_0\alpha_1\le4623.
\tag{6.3}
\]

The polynomial discriminant obeys

\[
 0<T^2-4N\le1159^2-4=1,343,277.
\tag{6.4}
\]

The field discriminant divides the discriminant of the order
\(\mathbf Z[\alpha]\).  This proves (0.3) and a finite field set.

### Theorem 6.1 (general integral right-pentagon field bound)

Every exact rank-six full-trace K3 carrier satisfying Section 1 and having
a compact right-angled prescribed eigen-nef chamber has ground-field
discriminant at most \(1,343,277\).  No relation between the two free
right-pentagon moduli is required.  The proof uses only the integral facet
reflections, their full five-equation Gram closure, and the two signed
places; completeness, saturation, and exact relative maximality are
retained but not needed to enlarge the arithmetic argument.

## 7. Sharper integral-coordinate corollary

Suppose additionally that \(A,C\in\mathcal O\).  This is stronger than the
automatic statement \(4A,4C\in\mathcal O\).  Then
\(r=A-1\) and \(s=C-1\) are algebraic integers with distinguished values
positive and conjugates in \((-1,0)\).  Their nonzero integral norms give

\[
                                  r_0>1,\qquad s_0>1.
\tag{7.1}
\]

Using \(4x_3\in\mathcal O\), equation (4.4), and (3.4),

\[
 {1\over16}<x_3={s+1\over rs-1}
 \quad\Longrightarrow\quad
 r<16+{17\over s}<33.
\tag{7.2}
\]

Symmetrically \(s<33\), so

\[
                          1<A_0,C_0<34.
\tag{7.3}
\]

For the integral generator \(A\),

\[
 3\le\operatorname {Tr}(A)\le34,
 \qquad
 1\le N(A)\le33,
\tag{7.4}
\]

and hence

\[
                         |\operatorname {Disc}(k)|\le34^2-4=1152.
\tag{7.5}
\]

Exact enumeration of the box (7.4), with positive nonsquare polynomial
discriminant and roots in the intervals (2.2)--(2.3), gives 528 possible
trace--norm rows and 336 distinct squarefree radicands.  This proves (0.4).
The archived 336-element list is only a candidate field set; no existence
or realization count is asserted.

## 8. Poison packet against division by four

The sharper corollary cannot be silently promoted from pair-reflection
integrality.  Let

\[
 k=\mathbf Q(\sqrt5),\qquad
 \mathcal O_k=\mathbf Z[\phi],\qquad \phi^2=\phi+1,
\]

and consider the formal cyclic packet

\[
\begin{aligned}
 x_0&={\phi-1\over2},&
 x_1&={\phi-1\over4},&
 x_2&={1\over2}+2\phi,\\
 x_3&=8+14\phi,&
 x_4&=2+4\phi.&
\end{aligned}
\tag{8.1}
\]

It satisfies all five recurrences (3.3).  At the distinguished embedding,
the coordinates are approximately

\[
              (0.309,0.155,3.736,30.652,8.472),
\]

while at the conjugate embedding they are approximately

\[
              (-0.809,-0.405,-0.736,-0.652,-0.472).
\]

Every \(4x_i\) and hence every pair cycle \(4(x_i+1)\) lies in
\(\mathcal O_k\).  The directed five-cycle can also be chosen in the order;
with the acute distinguished sign it is

\[
                              c_5=-408-664\phi.
\tag{8.2}
\]

Nevertheless \(x_0,x_1,x_2\notin\mathcal O_k\), and in particular
\((x_0)_0<1\).  Thus the conclusion (7.1), and consequently the sharp
\(33/1152/336\) argument, genuinely requires the extra integrality of
\(A,C\).  The general quarter-integral theorem survives: all five numbers
still lie in the interval (0.1).

This is a poison **Gram-coordinate packet**, not a claimed K3 surface,
integral Cartan module, saturated mirror arrangement, or reflection-maximal
carrier.  It certifies precisely the invalid division-by-four step and
nothing broader.

## 9. Scope

Theorem 6.1 settles the independent-moduli compact right-pentagon field
gate.  It does not:

- prove that every exact \(n=2\) chamber is a right pentagon;
- treat noncompact polygons or compact Coxeter polygons with a different
  incidence graph;
- construct any carrier satisfying the exact hypotheses;
- bound a canonical K3 correspondence intersection
  \(H\mathbin\cdot\tau_\theta(H)\) for arbitrary chambers; or
- derive a general height theorem from reflection maximality alone.

The surviving \(n=2\) problem has therefore moved: one must either reduce
an arbitrary exact maximal eigen-nef chamber to a bounded collection of
incidence recurrences of this kind, or bound its length-at-most-five cyclic
coordinates without right-pentagon closure.
