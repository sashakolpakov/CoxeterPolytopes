# An explicit D4^2+A1^2 order-four carrier for the H10/38 section

## Positive theorem

Put \(s=t^2\), \(u=s-1\), and

\[
\begin{aligned}
 P(s)&=-3+\alpha u+\beta u^2,\\
 Q(s)&= 2-\alpha u+\gamma u^2+\delta u^3.
\end{aligned}
\]

Let \(X_{\alpha,\beta,\gamma,\delta}\) be the minimal resolution of

\[
 y^2=x^3+t^2P(t^2)x+t^3Q(t^2),
\]

and let

\[
 \sigma(x,y,t)=(-x,iy,-t).
\]

On a nonempty Zariski-open parameter set the singular fibers are

\[
 2I_0^*+2I_2+8I_1.
\]

For a very general member,

\[
 NS(X)=U\oplus D_4(-1)^2\oplus A_1(-1)^2.
\]

The primitive non-symplectic automorphism \(\sigma\) has order four. It
fixes both \(D_4\) diagrams vertex by vertex and exchanges the two
\(A_1\) components. Its Picard action is therefore the exact action of
Brandhorst--Hofmann row 0.4.4.6. After the certified Kähler/Weyl chamber
transport, this gives an explicit algebro-geometric realization of
R12-D4x2A1x2-p67-H10-38: a finite-volume noncompact Coxeter polytope in
\(H^{10}\) with 38 facets.

This is a new equation-to-Coxeter realization in the project. The abstract
polytope was already certified, and no worldwide priority for the abstract
polytope or the ambient order-four K3 family is claimed.

## Inverse-designed fibers

The parity of the coefficients gives

\[
 \sigma^2(x,y,t)=(x,-y,t),\qquad
 \sigma^*(dt\wedge dx/y)=-i\,dt\wedge dx/y.
\]

Write

\[
 E(s)=4P(s)^3+27Q(s)^2.
\]

The normalization at \(s=1\) gives

\[
 P(1)=-3,\qquad Q(1)=2,\qquad Q'(1)=-P'(1),
\]

and hence

\[
 E(s)=(s-1)^2R(s),\qquad \deg R=4.
\]

The leading coefficient at the double root is

\[
 R(1)=-9(\alpha^2-12\beta-12\gamma).
\]

When it is nonzero, the fibers over \(t=1\) and \(t=-1\) are exactly
\(I_2\). Their singular cubics are

\[
\begin{array}{c|c}
 t=1 &(x-1)^2(x+2),\\
 t=-1 &(x+1)^2(x-2).
\end{array}
\]

Thus \(x\mapsto-x\) exchanges the two resolved \(A_1\) components.

The full discriminant is

\[
 \Delta(t)=-16t^6E(t^2).
\]

At \(t=0\), the vanishing orders are \((2,3,6)\) whenever
\(E(0)\ne0\), so the fiber is exactly \(I_0^*\).

At infinity use \(v=1/t\), \(X=v^4x\), and \(Y=v^6y\). The transformed
coefficients have leading terms

\[
 a_4=\beta v^2+\cdots,\qquad a_6=\delta v^3+\cdots,
\]

and leading discriminant coefficient

\[
 4\beta^3+27\delta^2.
\]

Away from its zero locus the vanishing orders are again \((2,3,6)\), so
the fiber at infinity is a second \(I_0^*\). If \(R\) is squarefree and
its roots avoid the additive fibers and the zeroes of \(P\), its four
remaining \(s\)-roots lift to eight \(I_1\) fibers. The Euler count is

\[
 2\cdot6+2\cdot2+8=24.
\]

The open conditions are simultaneously nonempty. For

\[
 (\alpha,\beta,\gamma,\delta)=(0,1,0,1)
\]

one has

\[
 P=u^2-3,\qquad Q=u^3+2,
\]

and

\[
 R=31u^4-36u^2+108u+108.
\]

The exact witness data are

\[
\begin{aligned}
 \operatorname{disc}(R)&=-1900537572096,\\
 E(0)&=-5,\\
 4\beta^3+27\delta^2&=31,\\
 \operatorname{Res}(R,P)&=42849.
\end{aligned}
\]

All are nonzero, proving that the intended fiber stratum is nonempty.

## Why neither D4 diagram is permuted

This local calculation is essential: knowing only the abstract fiber types
would leave open a nontrivial triality action on each \(D_4\).

At \(t=0\), introduce the Tate-root coordinate

\[
 \xi=x/t.
\]

This is the first blowup chart of the actual branch plane. Substitution in
the Weierstrass cubic gives the exact identity

\[
 x^3+t^2P(t^2)x+t^3Q(t^2)
 =t^3\bigl(\xi^3+P(t^2)\xi+Q(t^2)\bigr).
\]

The three outer components of the resolved \(I_0^*\) fiber are indexed by
the three roots of

\[
 \xi^3+P(0)\xi+Q(0)=0.
\]

On the stated open set its discriminant is nonzero. Moreover

\[
 \sigma(\xi)=\frac{-x}{-t}=\xi.
\]

Each root and hence each outer component is fixed individually. The central
component is the unique central vertex and is fixed as well.

At infinity, in the coordinates above, \(\sigma\) acts by

\[
 (X,Y,v)\longmapsto(-X,iY,-v).
\]

Now put \(\xi=X/v\). The outer components are indexed by the roots of

\[
 \xi^3+\beta\xi+\delta=0,
\]

whose discriminant is

\[
 -4\beta^3-27\delta^2\ne0.
\]

Again

\[
 \sigma(\xi)=\frac{-X}{-v}=\xi,
\]

The transformed Weierstrass equation has the exact form

\[
Y^2=X^3+v^2\widehat P(v^2)X+v^3\widehat Q(v^2),
\]

and substituting \(X=v\xi\) factors its branch cubic as

\[
v^3\bigl(\xi^3+\widehat P(v^2)\xi+\widehat Q(v^2)\bigr).
\]

Thus \((v,\xi)\mapsto(-v,\xi)\) is again the literal first blowup chart,
not merely a Tate-root heuristic. The three branches hit the exceptional
divisor at distinct, individually fixed points. The double-cover
resolution therefore fixes all three outer components and the unique
central component. Both \(D_4\) diagrams are fixed vertex by vertex;
neither triality group is used.

## Exact generic Néron--Severi lattice

The zero section, fiber, and reducible-fiber components give

\[
 L=U\oplus D_4(-1)^2\oplus A_1(-1)^2\subset NS(X).
\]

It has rank 12, signature \((1,11)\), determinant \(-64\), and Smith
factors

\[
 1^6\,2^6.
\]

The parameters are recovered exactly from the normalized coefficients:

\[
 \alpha=P'(1),\qquad
 \beta=\frac{P''(1)}2,\qquad
 \gamma=\frac{Q''(1)}2,\qquad
 \delta=\frac{Q'''(1)}6.
\]

Their Jacobian with respect to
\((\alpha,\beta,\gamma,\delta)\) is exactly \(I_4\).

The fiber-type divisors \(\{0,\infty\}\) for \(I_0^*\) and
\(\{1,-1\}\) for \(I_2\) fix the base up to the four transformations
\(t\mapsto ct\) and \(t\mapsto c/t\) with \(c^2=1\). A
short-Weierstrass homothety preserving \(P(1)=-3,Q(1)=2\) satisfies
\(\lambda^4=\lambda^6=1\), hence \(\lambda^2=1\). Thus the normalized
coordinate stabilizer is finite.

The marked fiber and section give integral classes \(F,O\in NS(X)\) with
\(F^2=0\), \(O^2=-2\), and \(F\cdot O=1\). They are locally constant in
a connected fiber of the forgetful map. Fixed \(F,O\) recover the pencil
\(|F|\) and the section, while \(H^0(X,T_X)=0\) makes the K3 automorphism
group discrete. After the finite coordinate stabilizer above, the
marked-fibration-to-K3 map therefore has zero-dimensional fibers. Hence
the displayed four parameters give four genuine K3 moduli.

For a primitive order-four K3 action the transcendental lattice is a
\(\mathbf Q(i)\)-module, so the Picard rank is even. A strict generic jump
from \(\rho\ge12\) would therefore give \(\rho\ge14\). The corresponding
order-four period ball has dimension at most three, contradicting this
four-dimensional family. Hence the very general Picard rank is 12.

It remains to rule out a proper saturation of \(L\). Shioda--Tate gives
Mordell--Weil rank zero. In the height formula, a section meeting an
allowed nonidentity component of an \(I_0^*\) fiber contributes at most
one, and one meeting the nonidentity component of an \(I_2\) contributes
one half. Therefore a nonzero torsion section would satisfy

\[
 0=h(P)=4+2(P\cdot O)-\sum_v\operatorname{contr}_v(P)
       \ge4-(1+1+1/2+1/2)=1,
\]

a contradiction. Mordell--Weil torsion is trivial, the trivial lattice is
saturated, and

\[
 NS(X)=L
\]

for a very general member.

The rational quadruple above certifies the nonempty fiber open; it is not
claimed to be an individually Picard-rank-12 rational point without an
additional arithmetic reduction argument.

## Exact Picard action and BH handoff

Let \(R_+\) and \(R_-\) denote the two \(I_2\) roots. The component
calculation gives identity on \(U+D_4(-1)^2\) and

\[
 R_+\longleftrightarrow R_-.
\]

Consequently

\[
\begin{aligned}
 NS(X)^\sigma&\cong U+D_4(-1)^2+\langle-4\rangle,\\
 (NS(X)^\sigma)^\perp&\cong\langle-4\rangle.
\end{aligned}
\]

The fixed lattice has rank 11, signature \((1,10)\), determinant 64, Smith
factors \(1^6\,2^4\,4\), and discriminant quadratic histogram

\[
 \{0:16,\ 1:16,\ 3/4:12,\ 7/4:20\}.
\]

The coinvariant generator \(R_+-R_-\) has square \(-4\), divisibility two,
and is rootless. The fixed-plus-coinvariant lattice has index two in
\(NS(X)\), with the literal glue

\[
 \frac{(R_++R_-)+(R_+-R_-)}2=R_+.
\]

These invariants, the pinned discriminant-form isometry, and all primitive
glue pairs match Brandhorst--Hofmann row 0.4.4.6. The existing integral
lifting certificate therefore identifies the full Néron--Severi action
pairs, not only their rational characters.

The full characteristic polynomial is

\[
 (\lambda-1)^{11}(\lambda+1)(\lambda^2+1)^5.
\]

The corresponding BH fixed-locus fingerprints are

\[
\begin{array}{c|c}
 \sigma &8\text{ isolated points}+2\text{ rational curves},\\
 \sigma^2&3\text{ rational curves}+1\text{ genus-2 curve}.
\end{array}
\]

The BH row is effective and supplies an invariant Kähler class. Weyl
transport takes its Kähler chamber to Roulleau's complete nef chamber, and
the exhaustive weighted chamber census identifies the transported action
with the archived p67 involution.

## Literal finite-volume conclusion

The hyperbolic statement uses the full Vinberg criterion, not signature
alone. The archived complete chamber was exhaustively checked at every
elliptic codimension-one edge and every maximal parabolic neighborhood.
There are no failed edge or parabolic audits:

| invariant | value |
|---|---:|
| hyperbolic dimension | 10 |
| facets | 38 |
| finite vertices | 526 |
| cusps | 141 |
| elliptic edges | 4,518 |
| failed Vinberg audits | 0 |

The presence of cusps makes the chamber noncompact; the complete edge and
parabolic census proves finite volume. Thus the equation above realizes a
genuine finite-volume noncompact hyperbolic Coxeter polytope.

The canonical digest of the complete finite-vertex, ideal-component, edge,
and parabolic audit is
79c800f246b68e6eff1908dd1922abdf7b31fa3ac436e43920fce8b22a47cf03.

## Lineage and reproduction

The parity construction lies inside
[Artebani--Sarti Example 6.3](https://arxiv.org/abs/1102.4436). The present
normal form is the inverse-designed specialization which inserts
\(I_0^*/D_4\) fibers at zero and infinity and an exchanged \(I_2/I_2\)
pair at \(t=\pm1\).

Artifacts:

- [certificate](../scripts/certify_ud4x2a1x2_order4_weierstrass_family.py)
- [machine archive](../data/enumerations/ud4x2a1x2_order4_weierstrass_family.json)
- [tests](../tests/test_ud4x2a1x2_order4_weierstrass_family.py)

Rebuild and test with:

    PYTHONPATH=.:src python3 scripts/certify_ud4x2a1x2_order4_weierstrass_family.py --check-archive
    PYTHONPATH=.:src pytest -q tests/test_ud4x2a1x2_order4_weierstrass_family.py

For an external replay, add
--brandhorst-hofmann-zip /path/to/K3Groups.zip to the certificate command.
