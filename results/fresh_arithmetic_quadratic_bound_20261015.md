# Exact triangle restrictions for the quadratic-five arithmetic field

The [checker](../scripts/fresh_arithmetic_quadratic_bound_20261015.py) and [certificate](fresh_arithmetic_quadratic_bound_20261015.json) give a complete necessary triangle test at an \(E_7\) or \(E_8\) component of a finite vertex in a compact arithmetic Coxeter polytope with ground field \(K=\mathbb Q(\sqrt5)\). Exactly ten marked pairs in each component can occur on a triangular two-face. These restrictions remain valid when the vertex has other orthogonal spherical components.

They also reject all 23 **specific** E7/E8 quadrilateral recipes in the frozen [geometric small-face barrier](fresh_compact29_small_faces_20261015.md): twenty violate the cyclic-product field, and three violate arithmetic integrality. This is not an exhaustive quadrilateral classification. No improved dimension bound follows from this package alone.

All earlier files remain unchanged. Generation and saved replay pass. The [independent review](fresh_arithmetic_quadratic_bound_20261015_review.md) also passes, using separate Fraction arithmetic for all 8,748 profiles, the spectator argument and explicit rescue, and all 23 quadrilateral field/integrality failures.

## Normalization and complete finite-contact alphabet

All roots have squared length two. The E-type numbering is the chain \(0-1-\cdots-(r-2)\), with an extra edge \(2-(r-1)\). Its Gram matrix is \(G\), and \(J=G^{-1}\) is rational and entrywise positive. Write \(b_i=-\langle x,e_i\rangle\ge0\) for a new root's contacts and \(q=b^T Jb\) for its projection cost.

For an arithmetic chamber, cyclic products of this norm-two Gram are algebraic integers in the ground field, and conjugation away from its physical embedding makes the Gram positive semidefinite. These are Vinberg's necessary conditions; the field criterion and determinant consequences are recalled in [Vinberg, §7.1](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt). The integrality normalization is the cyclic-product condition on \(2G\) when the source uses unit-length normals; see [Vinberg's criterion, Theorem 2.1](https://msp.org/agt/2025/25-4/agt-v25-n4-p14-s.pdf).

For a finite label \(m\), the two-cycle gives \(\mathbb Q(\cos(2\pi/m))\subseteq K\). The real cyclotomic degree formula leaves \(m=2,3,4,5,6,8,10,12\); the fields at 8 and 12 are respectively \(\mathbb Q(\sqrt2)\) and \(\mathbb Q(\sqrt3)\), so they are excluded. Thus the complete list is

\[
m\in\{2,3,4,5,6,10\},\qquad
b_i\in\{0,1,\sqrt2,\phi,\sqrt3,\sqrt{\phi+2}\},
\quad \phi=(1+\sqrt5)/2.
\]

There is no finite-label cutoff here. Any two nonzero contacts close a cycle through the unique old E-tree path, whose edge magnitudes are all one. Hence their product belongs to \(K\), so all nonzero contacts have a common square class. The classes of \(2,3,\phi+2\) are distinct and nontrivial: a positive rational is a square in \(K\) only if it is a rational square or five times one, while \(N(\phi+2)=5\) excludes every class containing that factor from being a square. The alternatives are therefore a mixed \(\{1,\phi\}\) profile, or a uniform nonzero magnitude \(\sqrt2,\sqrt3,\sqrt{\phi+2}\).

The inverse diagonal minimum is \(3/2\) for E7 and 2 for E8. For the three uniform alternatives, conjugate projection cost is at least this minimum times \(2,3,3-\phi\), respectively. Each is strictly greater than two, including \(\tfrac32(3-\phi)>2\). Thus they violate conjugate positivity. Only contacts \(0,1,\phi\) remain.

## Exact one-root census

For each of the \(3^7\) and \(3^8\) profiles the script evaluates \(q\) and \(q^\sigma\) in exact Fraction arithmetic in the basis \((1,\phi)\). A full E-plus-root star is Lorentzian precisely when \(q>2\); its conjugate is positive definite precisely when \(q^\sigma<2\). Its determinant cannot become zero under conjugation if it was nonzero physically.

| Component | Profiles tested | Lorentzian with positive conjugate | Profiles realizing at least one triangle | Covered marked pairs |
|---|---:|---:|---:|---:|
| E7 | 2,187 | 242 | 5 | 10 of 21 |
| E8 | 6,561 | 286 | 3 | 10 of 28 |

Deleting node \(a\) changes the projection cost to

\[
q_a=q-\frac{(Jb)_a^2}{J_{aa}}.
\]

A marked pair \(\{a,b\}\) has the other two finite triangle corners exactly when both deletion costs are less than two. The checker evaluates all of these inequalities without rounded eigenvalues. Together with the old spherical corner and the full Lorentzian star, they give an actual compact triangle in the two-dimensional residual plane. The accepted Gram entries lie in \(\mathbb Z[\phi]\), so these local witnesses satisfy the required integral cyclic-product and conjugate-sign conditions. They are separate local realizations, not one simultaneous compact chamber.

The ten marked pairs are:

| Component | Pairs |
|---|---|
| E7 | 01, 02, 12, 23, 24, 25, 26, 34, 35, 45 |
| E8 | 01, 02, 12, 27, 34, 35, 36, 45, 46, 56 |

Here 01 denotes \(\{0,1\}\), and similarly elsewhere. For E7 the five successful contact profiles are single \(\pi/5\) contacts at nodes 0, 4, 5 or 6, and the double profile \((0,3),(6,5)\). For E8 they are single \(\pi/5\) contacts at 0 or 6, and \((0,3),(7,5)\).

## Orthogonal spectator components do not add triangle pairs

The E7/E8 block need not be the whole finite vertex diagram. Let \(Z\) be the union of the other spherical components, all white for the two-face in question, and let \(c\ge0\) be the new root's projection cost onto \(Z\). This requires care: an isolated E7-plus-root star need not itself be Lorentzian even though the full star is.

If the E projection already has \(q>2\), the full finite-corner inequalities \(q_a+c<2\) imply \(q_a<2\), so its pair is in the recorded list. If all E contacts vanish, deleting a marked E node changes no projection cost, so the full star cannot be indefinite while either such corner is positive definite.

The exact census finds no other E8 case. For E7 there is exactly one remaining profile: a single \(\pi/3\) contact at leaf 5, with \(q=3/2\). In this case full Lorentzian signature and conjugate positivity require

\[
c>1/2,\qquad c^\sigma<1/2.
\]

Thus the spectator cost is nonrational. If all spectator components attached to the new root used only bonds 3 and 4, their Gram entries would be in \(\mathbb Q(\sqrt2)\), whereas each Schur cost belongs to the cyclic ground field \(K\), by principal determinant ratios. Their intersection is \(\mathbb Q\), a contradiction. Consequently a bond 5 occurs in an attached component of an actual finite corner retaining leaf 5. Bonds at least six cannot occur there: that connected component contains leaf 5, the new root, and at least one spectator, whereas such finite bonds occur only in isolated rank-two components. A connected finite component containing a bond 5 is H3 or H4, of rank at most four.

It follows that the remaining E7 component containing leaf 5 has at most two nodes. Deleting the nodes 0 through 6 leaves that component with sizes

\[
6,5,3,2,1,0,6.
\]

Hence each marked node lies in \(\{3,4,5\}\), including the option of deleting leaf 5 itself. All three such pairs already belong to the ten-pair list. This proves the ambient assertion without assuming the other vertex components are absent.

There is a genuine spectator rescue, so omitting this argument would be incorrect. Adjoin an A1 spectator and give it a \(\pi/5\) contact to the new root, in addition to the \(\pi/3\) contact at E7 node 5. Then

\[
q+c=\frac32+\frac{\phi^2}{2}>2,
\qquad (q+c)^\sigma<2.
\]

The exact finite-deletion set is \(\{3,4,5\}\), as predicted. The certificate includes this example and the component-size calculation.

## The old quadrilateral barrier is not an arithmetic barrier

The script reads, without modifying, the 23 literal E7/E8 quadrilateral recipes from the frozen geometric barrier. In twenty cases a fundamental cycle through the dotted edge lies outside \(K\), certified by exact algebraic field membership. This cycle uses a nonzero new-to-old contact, the unique old E-tree path, the dotted edge, and the new-to-new edge; extra chords do not invalidate it.

The remaining three are tree-type cases whose dotted two-cycles are:

| Old recipe | Dotted square | Primitive integer minimal polynomial |
|---|---|---|
| E7, marks 03 | \(14/3+2\sqrt5\) | \(9X^2-84X+16\) |
| E7, marks 56 | \(5/2+9\sqrt5/10\) | \(5X^2-25X+11\) |
| E8, marks 03 | \(9/2+2\sqrt5\) | \(4X^2-36X+1\) |

Their primitive irreducible polynomials are not monic, so none is an algebraic integer. They therefore violate a necessary condition for an arithmetic chamber even where conjugate positivity holds. All 23 literal witnesses are rejected; **other quadrilateral completions have not been excluded**.

## What remains unresolved

This supplies triangle incidences that any hypothetical quadratic-five arithmetic chamber must obey, and removes an inappropriate reuse of a nonarithmetic local counterexample bank. A stronger dimension bound would need a complete treatment of relevant quadrilateral extensions or a globally justified face-count argument using these restrictions. Neither is established here. In particular, “not a triangle” does not mean “not a compact two-face.”

```sh
python3 scripts/fresh_arithmetic_quadratic_bound_20261015.py
python3 scripts/fresh_arithmetic_quadratic_bound_20261015.py --verify
```

Both commands pass. The source, certificate and report are final and frozen; no earlier source or result was changed.
