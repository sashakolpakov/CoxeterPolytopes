# Independent hostile audit: universal weak Chow action and leverage no-go

## Verdict

**PASS, with the owner's marked weak-category and imported-bound scope.**

The frozen package correctly constructs a faithful \(W(P)\)-action on a
rank-\(R\) divisor Tate summand, returns the distinguished chamber exactly,
and proves only the intended leverage no-go. In particular, the apparently
reversed coefficient formula

\[
                     C_AJC_B=C_{BA}
\]

has the right orientation: it says that the correspondence for \(A\),
followed by the correspondence for \(B\), acts as \(BA\). Consequently
\(\Gamma_g\circ\Gamma_h=\Gamma_{gh}\), so the construction is an action, not
an anti-action. The independent replay uses a noncommuting pair; reversing
the order fails.

No owner byte was edited. The terminal owner SHA-256 pins are

~~~text
df38b26e0cffd8a391834be922869153af5badf41a632e74a48b472b332e101a  report
d3a19cfe7069b811721a77a624952d92ab5bc7c5574d9e332d2ff8c7f0205eba  certifier
d75897bbf8d80a662a7b9d3aa3ae729052494e20633d835566db70870bb0fc9e  archive
6da7b9e4590d7e531a7f3502ba2c971e806cbd3f259cb4722169553b7532c802  tests
~~~

The separate replay is in the [audit script](../scripts/audit_universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile.py),
[audit archive](../data/enumerations/universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile_audit.json),
and [audit tests](../tests/test_universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile_audit.py).

## 1. The stable trace lattice is legitimate

Let \(\Gamma=W(P)\) and start with an \(\mathcal O_k\)-lattice \(L\) whose
stabilizer meets \(\Gamma\) in finite index. Take the normal core

\[
 K=\bigcap_{\gamma\in\Gamma}
       \gamma\bigl(\Gamma\cap\operatorname {Stab}(L)\bigr)\gamma^{-1}.
\]

Then \(K\triangleleft\Gamma\) has finite index and preserves \(L\). If \(T\)
is a finite set of representatives for \(\Gamma/K\), put

\[
                         M=\sum_{t\in T}tL.
\]

For \(\gamma\in\Gamma\), write \(\gamma t=t'k\), with \(t'\in T\) and
\(k\in K\). Hence \(\gamma(tL)=t'L\): \(\Gamma\) permutes the finitely many
summands and preserves \(M\). Because every \(tL\) is a full lattice, so is
\(M\). This closes the finite-index-to-full-action step.

The values \(q(M,M)\) lie in a fractional ideal. Its inverse, multiplied by
the codifferent, contains a positive rational number, hence a totally
positive \(\beta\) for which

\[
             B(x,y)=\operatorname {Tr}_{k/\mathbf Q}(\beta q(x,y))
\]

is integral. Positive rational scaling also puts each of the finitely many
oriented root representatives in \(M\); it changes neither its reflection nor
its \(\sigma _0\)-oriented wall. At the real places, \(B\) has one positive
direction from \(\sigma _0\) and only negative directions elsewhere, so its
signature is \((1,R-1)\).

The executable audit includes a rational involution which does not preserve
\(\mathbf Z^2\), constructs \(L+gL\), and verifies in an integral basis that
the summed lattice is stable and its rescaled Lorentzian Gram matrix is
integral.

## 2. The \(4R-3\) embedding count is exact

After rational congruence diagonalization, write

\[
 G\sim_{\mathbf Q}\langle a_0,-a_1,\ldots,-a_{R-1}\rangle,
 \qquad a_i>0.
\]

The identity

\[
 a_0=\left({a_0+1\over2}\right)^2-
     \left({a_0-1\over2}\right)^2
\]

uses the unique positive ambient coordinate and one negative coordinate.
Each \(a_i\), \(i>0\), is a sum of four rational squares and uses four new
negative coordinates. Clearing one common denominator therefore gives

\[
 j^tI_{1,b}j=D^2G,\qquad
 b\le 1+4(R-1)=4R-3.
\]

The audit independently realizes
\(\langle7/3,-5/6,-11/10\rangle\) in \(I_{1,9}\), clears denominator \(30\),
and verifies the scaled Gram identity entry by entry. This tests both the
positive difference-of-squares block and two disjoint four-square blocks;
it does not call the owner embedding routine.

## 3. Literal Chow composition, transpose, and action orientation

Let \(e_i\) be the divisor basis, \(J=((e_i,e_j))\), and use column vectors
for endomorphisms. For

\[
 Z_A=\sum_i e_i^*\times A(e_i),
\]

the coefficient matrix in the \(e_i\times e_j\) basis is

\[
                         C_A=J^{-1}A^t.
\]

Push-pull on a divisor column \(x\) gives

\[
                         C_A^tJx=Ax.
\]

If \(Z_A\) is followed by \(Z_B\), intersection of the two middle divisor
factors gives

\[
 C_AJC_B=J^{-1}A^tB^t=J^{-1}(BA)^t=C_{BA}.             \tag{3.1}
\]

Thus, under the standard convention that \(\Gamma_g\circ\Gamma_h\) means
first \(h\), then \(g\), equation (3.1) yields

\[
                         \Gamma_g\circ\Gamma_h=\Gamma_{gh}.
\]

There is no orientation defect.

For the orthogonal projector

\[
 \Pi=j(j^tJj)^{-1}j^tJ
\]

and the zero-complement lift

\[
 L_g=jg(j^tJj)^{-1}j^tJ,
\]

one has \(L_gj=jg\), \(L_gL_h=L_{gh}\), and \(L_1=\Pi\). A native
reflection is self-adjoint for the source form, so \(L_s^tJ=JL_s\). This is
equivalent to \(C_{L_s}=C_{L_s}^t\), proving self-transpose at the literal
cycle level, while \(L_s^2=\Pi\) proves it is an involution of the summand.

The independent \((2,3,\infty)\) replay uses an explicit rank-three injection
into \(I_{1,5}\). It verifies the \(m=2\) and \(m=3\) relations, the
nontrivial parabolic word, self-transpose, and both orders of a noncommuting
pair. Only the correct order satisfies (3.1).

## 4. The distinguished eigenslice returns the chamber exactly

Restriction of scalars gives the orthogonal real decomposition

\[
 (M_{\mathbf Z}\otimes\mathbf R,B)
   =\bigoplus_{\sigma:k\hookrightarrow\mathbf R}
     (M_\sigma,\sigma(\beta)q_\sigma).
\]

Multiplication by \(k\) is \(B\)-self-adjoint and splits this space into its
real embedding summands. If \(x\in E_{\sigma _0}\), all other components of
\(x\) vanish, so for every transported wall \(D_i=j(r_i)\),

\[
 (x,D_i)_X=D^2\sigma _0(\beta)q_{\sigma _0}(x,r_i).    \tag{4.1}
\]

The factor is positive, common to all walls, and also scales the metric on
the slice. Therefore (4.1) preserves the chosen half-spaces, normalized
angles, ultraparallel distances, incidences, and the selected positive-cone
component.

This step is independently replayed over
\(k=\mathbf Q(\sqrt2)\) with \(q=\langle\sqrt2,-1\rangle\). Its trace form
has signature \((1,3)\); multiplication by \(\sqrt2\) has a Lorentzian
\(+\sqrt2\) eigenspace and a negative-definite \(-\sqrt2\) eigenspace. Two
oriented walls have exact trace/slice pairings \(1\) and
\(2+2\sqrt2\), both positive. This is a genuine nontrivial field-slice test,
not the owner's rational triangle calibration. The replay also uses the
field-linear reflection which fixes the first field coordinate and negates
the second; it is a trace isometry, an involution, and commutes exactly with
multiplication by \(\sqrt2\).

## 5. The number \(2997\) is imported, not derived

The construction gives \(b\le4R-3\) with \(R=d(n+1)\). The independently
pinned arithmetic-bound package supplies the separate published inputs

\[
                  n+1\le30,\qquad d=[k:\mathbf Q]\le25.
\]

Only after those inputs are substituted does one obtain

\[
 R\le750,\qquad b\le4\cdot750-3=2997.
\]

The owner labels this dependency and circularity correctly. The bound is
valid for finite-volume arithmetic inputs and does not become an
algebraic-geometric proof of either classical inequality.

## 6. The leverage no-go survives hostile pressure

Internally, the carrier supplies only

\[
                 d(n+1)=R\le\rho(X)=b+1.
\]

Because \(b\) is deliberately allowed to grow with \(R\), this is not a
ceiling. On the generic construction \(b+1=4R-2\), so the inequality is
tautological for every \(R\).

Moreover, divisor-product cycles realize all of
\(\operatorname {Mat}_{\rho(X)}(\mathbf Q)\). Already on a fixed rank-two
Tate part, the matrices

\[
 A_D=\begin{pmatrix}0&D\\1&0\end{pmatrix},\qquad A_D^2=DI,
\]

give actions of the infinitely many real quadratic fields
\(\mathbf Q(\sqrt D)\) as \(D>1\) ranges over squarefree integers. Hence an
abstract field action has no field-discriminant or coefficient-height
scarcity. The owner's particular trace-compatible self-adjoint action may
vary with the input; the point is that the weak Chow category itself imposes
no independent bounded complexity.

This no-go does not deny the classical arithmetic bounds, and it does not
exclude leverage from finite volume, maximality, volume, effective/nef/MBM
geometry, or a canonical bounded carrier. It says precisely that weak
Chow-carrier existence alone cannot provide that leverage.

## Reproduction and presentation checks

~~~bash
PYTHONPATH=src:. python3 \
  scripts/audit_universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile.py \
  --check

PYTHONPATH=src:. pytest -q \
  tests/test_universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile_audit.py \
  tests/test_universal_weak_chow_action_carrier.py \
  tests/test_unconditional_maximal_arithmetic_reflection_chamber_finiteness.py
~~~

Both the owner and audit Markdown have balanced display delimiters, no
forbidden control characters, and convert successfully through Pandoc.
