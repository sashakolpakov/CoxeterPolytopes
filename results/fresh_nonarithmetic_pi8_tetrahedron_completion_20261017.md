# A compact tetrahedral face in the 58-wall H9 chamber

**Exact generation and fresh saved replay PASS.** Three compatible walls extend the frozen 55-wall chamber to a coherent 58-wall Coxeter chamber. Its actual face supported by
\[
S'=[0,1,2,3,51,55]\cong H_4+I_2(8)
\]
is exactly a compact Coxeter tetrahedron. This is a local construction result; neither finite volume of the ambient chamber nor closure of the original \(H_4+A_2\) face is claimed.

The certificate is [the executable source](../scripts/fresh_nonarithmetic_pi8_tetrahedron_completion_20261017.py) and [its saved data](fresh_nonarithmetic_pi8_tetrahedron_completion_20261017.json). It binds the frozen [complete endpoint locus](fresh_nonarithmetic_h3_edge730_domain_20261017.md), the [complete pi/5 tetrahedron catalogue](fresh_h9_pi5_tetrahedron_templates_20261017.md), and the prior 55-wall geometry. The implementation uses its frozen outward rational interval arithmetic, with exact equalities supplied by the displayed projection and quadratic recipes. An independent reconstruction is being performed separately.

## The three walls

Wall55 is the first admissible endpoint \(r_8\) from the complete old-edge \([7,30]\) locus. It is right to old roots \(0,1,2,3,5,6,7,30\), has angle \(\pi/8\) with51, and is strictly dotted to every other old wall.

Let \(u=56\), \(v=57\). Both are right to every root in \(S'\). Their remaining finite contacts are

| New wall | Additional right contacts | Additional nonright contacts |
|---|---|---|
| \(u\) | old6, old7 | old30 at \(\pi/3\) |
| \(v\) | old30 | old7 at \(\pi/3\), \(u\) at \(\pi/4\) |

Every other new–old or new–new contact is strictly dotted. The largest dotted products, with all roots normalized to square2, are approximately −2.03298 for55, −3.84841 for56, and −2.10695 for57. These decimals only summarize the strict rational interval certificates.

The four perpendicular sides of the new face are \([7,30,56,57]\), with Gram matrix
\[
\begin{pmatrix}
2&-\varphi&0&-1\\
-\varphi&2&-1&0\\
0&-1&2&-\sqrt2\\
-1&0&-\sqrt2&2
\end{pmatrix},\qquad \varphi=(1+\sqrt5)/2.
\]
This is literal template9 of the frozen catalogue, with labels \([5,2,3,3,2,4]\) in the catalogue's pair order. The search that suggested it does not enter the proof.

## Exact realization and common interior

For numerical conditioning, the construction works in the ten-root basis
\([0,1,2,3,51,r_8,7,30,5,6]\). Its literal Gram matrix, the physical normal coordinates, and all conversion data needed to reconstruct the roots are saved.

Project the prescribed contact rows \((0,-1)\) and \((-1,0)\) onto the positive plane spanned by old7,30, obtaining \(p_u,p_v\). Put \(r_u=2-p_u^2\), \(r_v=2-p_v^2\), \(d=-\sqrt2-p_u\cdot p_v\), and \(\gamma^2=d^2/r_u-r_v\). Let \(k_0^2=-1,k_1^2=1,k_0\cdot k_1=0\) span the complement of the eight fixed roots, and set
\[
A=(k_0+k_1)/\sqrt2,\quad B=(-k_0+k_1)/\sqrt2.
\]
Then
\[
u=p_u+\sqrt{r_u/2}(tA+B/t),
\]
\[
v=p_v+(d/\sqrt{2r_u}+\gamma/\sqrt2)tA
 +(d/\sqrt{2r_u}-\gamma/\sqrt2)B/t.
\]
The parameter is the unique positive root of \(u\cdot n_6=0\), approximately3.8645469727. The saved quadratic has negative leading coefficient and positive constant coefficient. These identities give both norms2 and all prescribed products exactly; interval checks certify every remaining strict sign.

All three additions are strict at the same persistent old finite vertex \(V\). The explicit point \(X=V+H\), where \(H\) is the old strict center, is future timelike and strict against all58 walls. For every wall \(n_i\), the certificate checks
\[
X_i=X-\tfrac12(X\cdot n_i)n_i.
\]
It is future timelike, right to \(n_i\), and strict against every other wall. Thus each of the58 inequalities defines an actual facet, rather than a redundant halfspace.

## Entire-face closure

The four dual tetrahedron vertices are computed from the inverse of the displayed four-side Gram. Each is future timelike and satisfies all58 inequalities. Its incident set consists of all six roots of \(S'\) and exactly three tetrahedron sides; all four rank9 incident Grams are positive definite. Every nonincident inequality is strict.

The four-side simplex in \(S'^\perp\) is compact. Since each of its four vertices satisfies every remaining ambient inequality, their convex hull does too. Conversely, the actual \(S'\) face is contained in this four-side simplex. Hence the whole actual face equals the compact tetrahedron; this conclusion does not merely count four isolated local corners.

The original support \(S_6=[0,1,2,3,5,6]\) defines a different face. Both56 and57 are right to its H4 factor, and their A2 projections have residual square
\[
2-\tfrac23(g_5^2+g_5g_6+g_6^2)<0.
\]
For56, \(g_6=0\) and \(g_5<-2\); for57 both products are below−2. Their negative products with the future finite point \(V\in S_6^\perp\) fix the future orientation. Therefore both residual projections are future timelike and are strict on the entire future part of the original face. That face changes only through wall55. Closing it, and closing the remaining ambient exterior regions, requires further work.

Generation73679 and fresh saved replay34858 both exited0. The source and data are frozen after syntax, direct source-hash, and report-link checks. No complete ambient cone enumeration is included in this package.
