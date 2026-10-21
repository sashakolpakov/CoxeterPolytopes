# Independent hostile audit: the compact (2,4,6) exact Property-4 counterexample

**Verdict (2026-08-11): PASS at the owner's exact scope.** The triangle is
an arithmetic compact Coxeter chamber over \(\mathbf Q\), its side group has
no proper discrete reflection overgroup, and all three currently declared
**exact-return** carrier species are empty. This is a counterexample to
essential surjectivity of that three-species full-nef compiler. It is not a
counterexample to the conditional dimension theorem after a carrier has been
supplied, and it does not rule out a new carrier species, a proper
invariant/eigen-nef slice for a rational-field input, or the bare
carrier-relative third row without exact chamber return.

This audit is independent of the owner executable: it imports neither the
owner module nor SymPy. It uses exact Fraction arithmetic, explicit
three-by-three matrix multiplication, a different index-bounded maximality
census, and a requirement-by-requirement reconstruction of the four-property
scope. It did not edit any owner byte.

## 1. Frozen owner bytes and raw reconstruction

The owner package was unchanged throughout this audit:

~~~text
860d82c50ecfe46bd12c0f2fba851e85c6d08f6bf5e6b9b6ef575c75f31006b0  results/compact_246_raw_maximal_property4_counterexample.md
c0afe27bacc218b1128ece87de8b7fa1f421a574f8b9039afa54947c2d73a64c  scripts/certify_compact_246_property4_counterexample.py
74fbb33d09a31af0a0af3ec1a65694732362ccc8e6eb6d799a4906a4debf2cae  data/enumerations/compact_246_property4_counterexample.json
089410474b73f09c134a10b00c33edc96c6ab30b9e1e6171b26be7eaaf3c34e0  tests/test_compact_246_property4_counterexample.py
~~~

Starting only from

\[
G=
\begin{pmatrix}
-12&6&0\\
6&-6&3\\
0&3&-2
\end{pmatrix},
\]

direct expansion gives

\[
\det G=36,\qquad
\text{LDL pivots }(-12,-3,1).
\]

Thus the signature is \((1,2)\). For the pairs
\((1,2),(2,3),(1,3)\), the normalized squared pairings are

\[
\frac12,\qquad\frac34,\qquad0,
\]

so the Coxeter labels are \(4,6,2\). The corresponding rank-two
determinants are \(36,3,24\), hence all vertices are finite. The chamber is
compact and

\[
\operatorname{area}(P)
=\pi\left(1-\frac12-\frac14-\frac16\right)
=\frac{\pi}{12}.
\]

The four nonzero Cartan coefficients are

\[
C_{12}=-1,\qquad C_{21}=-2,\qquad
C_{23}=-1,\qquad C_{32}=-3.
\]

The audit independently regenerates the three integral involutions

\[
s_1=\begin{pmatrix}-1&1&0\\0&1&0\\0&0&1\end{pmatrix},
\quad
s_2=\begin{pmatrix}1&0&0\\2&-1&1\\0&0&1\end{pmatrix},
\quad
s_3=\begin{pmatrix}1&0&0\\0&1&0\\0&3&-1\end{pmatrix}
\]

and checks \(s_i^tGs_i=G\) exactly. The cyclic Cartan products are \(2\)
and \(3\), so the minimal Vinberg field is \(\mathbf Q\). Geometrically,
the compact Coxeter triangle already makes \(W(P)\) a lattice; since it is a
subgroup of the rational arithmetic lattice \(O(G,\mathbf Z)\), the index is
finite. This proves the advertised arithmeticity without a classification
input.

## 2. The rational form is genuinely anisotropic

For the rationally similar form \(q=-G/2\), the quadratic polynomial is

\[
q(x,y,z)=6x^2-6xy+3y^2-3yz+z^2.
\]

If a nonzero rational zero existed, clearing denominators and dividing by
the common gcd would give a primitive integral zero. Modulo \(3\), the
equation first gives \(z^2=0\), hence \(3\mid z\). Write \(z=3z_1\), divide
the equation by \(3\), and reduce again. The remaining binary equation is

\[
2x^2+xy+y^2=0\pmod3.
\]

The independent nine-pair census has only \((x,y)=(0,0)\). Therefore
\(3\mid x,y,z\), contradicting primitivity. The form is anisotropic over
\(\mathbf Q\). This fact will obstruct a rank-three rational-surface host,
independently of the nef-cone argument.

## 3. Absolute reflection maximality

Suppose \(W(P)\) had a proper discrete reflection overgroup \(W(Q)\) of
index \(N\ge2\). Mirrors of the overgroup refine the compact chamber, so
\(P\) is tiled by \(N\) compact \(Q\)-chambers and

\[
\operatorname{area}(Q)=\frac{\pi}{12N}\le\frac{\pi}{24}.
\]

A compact Coxeter polygon with four or more sides has area at least
\(\pi/6\): in the four-sided case at least one of the four labels exceeds
\(2\), and in five or more sides the all-right lower bound is already
\(\pi/2\). Hence \(Q\) is a triangle.

The least compact Coxeter-triangle area is \(\pi/42\), attained by
\((2,3,7)\). Therefore

\[
\frac1{12N}\ge\frac1{42},
\]

so only \(N=2,3\) require inspection. For ordered
\(2\le a\le b\le c\), exact area equality becomes

\[
\frac1a+\frac1b+\frac1c=1-\frac1{12N}.                 \tag{3.1}
\]

Equation (3.1) forces \(a\le3\) and then
\(b\le2/(1-1/(12N)-1/a)\); \(c\) is determined. Thus the executable loops
have derived finite bounds rather than a search cutoff. The complete census
is

\[
\begin{array}{c|c}
N&Q\\ \hline
2&(2,3,8)\\
3&\varnothing.
\end{array}
\]

For index two, connectedness makes \(P\) the reflection double of two
adjacent \((2,3,8)\) chambers. The three possible shared sides give

\[
\begin{array}{c|c}
\text{shared endpoint orders}&\text{union}\\ \hline
(2,8)&(3,3,4)\\
(2,3)&(\pi/8,\pi/8,2\pi/3),\ \text{not Coxeter}\\
(3,8)&\text{a quadrilateral}.
\end{array}
\]

None is \((2,4,6)\). Consequently \(W(P)\) is maximal among all discrete
hyperbolic reflection groups, which is stronger than either declared-host
or arithmetic-reflection maximality.

## 4. Both rational exact full-nef species are empty

Because the minimal field is \(\mathbf Q\), Chamber Fidelity in the frozen
contract requires literal equality

\[
P=\mathbf P(\operatorname{Nef}S),
\]

not a proper eigenslice. A full nef chamber of hyperbolic dimension \(2\)
forces \(\rho(S)=3\).

There are two independent contradictions.

First, a rank-three smooth rational surface is nonminimal. Its numerical
lattice is the odd unimodular lattice

\[
I_{1,2}=\langle1\rangle\oplus\langle-1\rangle^2,
\]

which is isotropic on \((1,1,0)\). A full-rank rational similarity from the
anisotropic trace space would preserve the existence of nonzero isotropic
vectors. Thus the required rational HOST cannot exist.

Second, every smooth projective rational surface with Picard rank at least
two has a nonzero nef isotropic class. If a rational MMP ends on a ruled
surface, pull back a fibre. If it ends on \(\mathbf P^2\), stop immediately
after the final inverse blowup at \(\mathbb F_1\) and pull back the class
\(H-E\). Pullback preserves nefness and square zero. Its ray is an ideal
point of the complete projectivized nef cone, whereas the \((2,4,6)\)
triangle has three finite vertices.

Therefore both current rational output spaces are empty:

\[
\operatorname{Car}_{\mathrm{direct}}(P)
=\operatorname{Car}_{\mathrm{weak\mbox{-}log}}(P)
=\varnothing.
\]

The weak-log decoration and semiampleness machinery never get a chance to
repair this: the obstruction is to the smooth rational host and literal
full-nef return themselves.

## 5. The exact-return minimalizable full-trace species is empty

Here the word **exact-return** does essential work. The third row of the
conditional transfer theorem by itself asks for a canonical-orthogonal or
finite-orbit-minimalizable trace block; it does not assert that its carrier
returns \(P\). After combining that row with Chamber Fidelity for the
rational field, literal full-nef equality again forces \(\rho(S)=3\).

The trace block also has rank three. A nontrivial contraction from this
surface has endpoint Picard rank at most two, so the entire rank-three block
cannot descend injectively. The exact carrier must already be minimal. On
the minimal branch the trace block lies in \(K_S^\perp\); because it spans
\(\operatorname{NS}(S)_{\mathbf Q}\), this forces

\[
K_S\equiv0.
\]

Let \(C\) be an irreducible physical negative curve defining a returned nef
facet. Adjunction gives

\[
C^2=2p_a(C)-2<0.
\]

Hence \(p_a(C)=0\), \(C^2=-2\), and \(C\) is a smooth rational curve. For
two distinct physical facets \(C,D\), the intersection \(C.D\) is a
nonnegative integer. A finite Coxeter angle requires

\[
\cos(\pi/m)=\frac{C.D}{2}<1,
\]

so \(C.D\in\{0,1\}\) and \(m\in\{2,3\}\). The target has labels \(4\) and
\(6\). This contradiction proves

\[
\operatorname{Car}^{\mathrm{exact}}_{\mathrm{min\mbox{-}trace}}(P)
=\varnothing.
\]

It does not prove that a bare rank-three trace block cannot occur on some
carrier when no exact chamber return is requested.

## 6. Exact four-property verdict

The result can be read directly against the four governing requirements:

| requirement | audit result on the current three species |
|---|---|
| 1. full trace/order HOST | the two rational-surface hosts fail by anisotropy versus \(I_{1,2}\) |
| 2. exact chamber | literal \(k=\mathbf Q\) full-nef return forces rank three; every rational full-nef cone has an ideal ray |
| 3. canonical/log/MMP fidelity | nontrivial rank-three descent is impossible; the minimal \(K\)-trivial lane permits only labels \(2,3\) |
| 4. essential surjectivity | false for the union of the direct rational, certified weak-log rational, and exact-return minimalizable full-trace species |

Thus the owner conclusion is exactly justified:

\[
\boxed{\text{raw absolute reflection maximality does not force a nonempty
current exact three-branch carrier list}.}
\]

The quantifiers must not be enlarged. The audit does not show that:

- the conditional AG dimension bounds are false;
- every conceivable algebraic-geometric carrier is impossible;
- a proper invariant/eigen-nef slice for \(k=\mathbf Q\) is impossible;
- a projected return in a larger Picard lattice is impossible;
- the bare carrier-relative third row is empty;
- a higher-dimensional variety or positive-characteristic carrier is
  impossible; or
- the abstract \((2,4,6)\) Coxeter chamber fails to exist.

Accordingly, this closes a no-go for the presently advertised exact-return
compiler, not for every possible future POLYTOPE-to-SURFACE functor.

## 7. Presentation and replay

The owner report is self-contained, distinguishes absolute reflection
maximality from host-relative maximality, explicitly marks the exact-return
strengthening of the third row, and gives the necessary escape clauses. Its
display mathematics converts successfully through Pandoc. The broad
common-class discussion in its Section 6 is correctly not promoted beyond
the frozen three-species theorem.

The independent artifacts are:

- [auditor](../scripts/audit_compact_246_property4_counterexample_independent_hostile.py);
- [exact archive](../data/enumerations/compact_246_property4_counterexample_independent_hostile_audit.json); and
- [regression tests](../tests/test_compact_246_property4_counterexample_independent_hostile_audit.py).

Run:

~~~text
python3 -m scripts.audit_compact_246_property4_counterexample_independent_hostile --check-output
pytest -q tests/test_compact_246_property4_counterexample_independent_hostile_audit.py
~~~
