# Explicit separation at classical irreducible cusps of a narrow facet

**Conditional geometric bounds, with all rank parameters explicit.** Let (P\subset\mathbb H^n) be a finite-sided finite-volume Coxeter polytope. Suppose a Nikulin narrow facet (F) contains an ideal vertex whose **full** cusp diagram is one irreducible affine classical diagram. Then the following bounds hold:

| Full cusp diagram | Bound on ambient dimension (n) |
|---|---:|
| (\widetilde A_{n-1}) | (n\le135) |
| (\widetilde B_{n-1}) | (n\le51) |
| (\widetilde C_{n-1}) | (n\le47) |
| (\widetilde D_{n-1}) | (n\le55) |

The proof supplies explicit parts of the missing narrow-facet argument. It does **not** prove a new unrestricted dimension bound, the full historical assertion (n<300), or a bound for every simple finite-volume polytope. Products of affine components and a chosen narrow facet having only finite vertices remain untreated. No novelty or strongest-known-bound claim is made.

The [exact checker](../scripts/fresh_nikulin_narrow_20261015.py) and [certificate](fresh_nikulin_narrow_20261015.json) verify closed inverse formulas, Lorentz-coordinate identities, energy inequalities and the boundary indexing with `Fraction` arithmetic. The proofs below, rather than a finite rank cutoff, establish the all-rank statements.

## Published input and geometric conventions

[Nikulin's ICM report](https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf), printed p.655, gives a facet (F) such that the hyperplanes of any two facets whose hyperplanes meet (H(F)) have distance less than (\operatorname{arccosh}7). Printed pp.657–659 describe terminating walls on cusp edges and illustrate affine (A), obtaining (n<140). We use the narrow-facet theorem and reconstruct the separation calculation independently. The [frozen source audit](dimension_bound_prokhorov_20261015_proof_review.md) records the exact primary provenance and the incomplete status of the full 300 reconstruction.

Use Lorentz signature ((n,1)), outward facet normals of square 2, and chamber inequalities ((h,\alpha)<0). Acuteness gives ((\alpha,\beta)\le0) for different walls. Every nonzero negative pairing has magnitude at least 1: intersecting walls give (2\cos(\pi/m)), (m\ge3), while parallel or disjoint walls give at least 2. Thus no finite-label cutoff is needed.

At the assumed cusp there are (n) incident walls, with normals (\alpha_0,\ldots,\alpha_{n-1}). Deleting node (i) gives a positive-definite rank-((n-1)) edge support (S_i). Its geometric edge has the original cusp as one endpoint. Choose any additional facet (x_i) at its other endpoint. Finite volume ensures that endpoint is finite or ideal, not beyond the hyperbolic closure. Consequently the Gram matrix on (S_i\cup\{x_i\}) is positive semidefinite. This remains true at a nonsimple ideal endpoint; no simplicity assumption is imposed there.

## 1. Schur localization, without endpoint-diagram enumeration

Write (J=G(S_i)^{-1}) and (p_j=-(x_i,\alpha_j)\ge0). Positive semidefiniteness gives

\[
p_{S_i}^{T}Jp_{S_i}\le2.
\]

All entries of the inverses below are nonnegative. Thus (J_{jj}>2) forces (p_j=0), because a nonzero Coxeter contact has (p_j\ge1).

The needed inverse formulas are elementary. For (A_r), in path order (a,b=0,\ldots,r-1),

\[
J_{ab}=\frac{(\min(a,b)+1)(r-\max(a,b))}{r+1}.
\]

For (B_r), number the short end 0 and the remaining path (1,\ldots,r-1), normalizing **all** roots to square 2:

\[
J_{00}=r/2,\quad J_{0j}=(r-j)/\sqrt2,\quad
J_{ij}=r-\max(i,j)\quad(i,j\ge1).
\]

For (D_r), use fork leaves 0,1, attached to 2, followed by the path to (r-1):

\[
J_{00}=J_{11}=r/4,\quad J_{01}=(r-2)/4,
\quad J_{0j}=J_{1j}=(r-j)/2\ (j\ge2),
\quad J_{ij}=r-\max(i,j)\ (i,j\ge2).
\]

Multiplication by their sparse Gram matrices proves these identities for every rank; all entries are positive. In (A_{n-1}), every diagonal except those within two positions of an end exceeds 2 when (n\ge10). In (B_r), the same assertion holds for (r\ge5); in (D_r), it holds for (r\ge9). The two allowed end diagonals in (B,D) are 2 and 1.

It follows that, when the retained components meet those rank thresholds, the full contact row of (x_i) is supported only on

\[
\{i-2,i-1,i,i+1,i+2\}.
\]

For affine (A) the indices are cyclic. For the other families they are interior path indices. The contact at the deleted node (i) itself is unrestricted except by acuteness and the Coxeter gap; it has **not** been assumed finite.

## 2. A Lorentz row identity and its uniform lower bound

Let (c=\sum\delta_j\alpha_j) be a future null cusp vector, with positive affine coefficients. Put

\[
p_j=-(x,\alpha_j),\quad q_j=-(y,\alpha_j),\quad
\sigma=\sum\delta_jp_j=-(x,c),\quad
\tau=\sum\delta_jq_j=-(y,c).
\]

Both (\sigma,\tau>0). A zero contact row would make a norm-two vector orthogonal to the full cusp span, whose orthogonal complement is the null line (\mathbb Rc), an impossibility. Let (E) be the squared Euclidean norm of the image of (x/\sigma-y/\tau) in (c^\perp/\mathbb Rc). Directly expanding its square gives the exact identity

\[
\boxed{-(x,y)=\frac{\sigma\tau}{2}E-\frac{\sigma}{\tau}-\frac{\tau}{\sigma}.}
\]

If all possible contact positions have the same affine coefficient (d), then the Coxeter gap gives (\sigma,\tau\ge d). When (E\ge4/d^2), the identity implies

\[
-(x,y)\ge\frac{d^2E}{2}-2.
\]

Indeed, (1/\sigma^2+1/\tau^2\le2/d^2) and (\sigma\tau\ge d^2). This removes any need for an upper bound on the unknown deleted-node contacts, or a limiting argument involving those contacts.

## 3. The path energy for affine (B,C,D)

Write (m=n-1). In a null basis (c,d_0) with ((c,d_0)=-1), together with orthonormal spatial coordinates (e_1,\ldots,e_m), one may use:

* (\widetilde C_m): (\alpha_0=\sqrt2(c-e_1)), (\alpha_j=e_j-e_{j+1}) for (1\le j<m), (\alpha_m=\sqrt2 e_m). The affine coefficients are (1/\sqrt2,1,\ldots,1,1/\sqrt2).
* (\widetilde B_m): (\alpha_0=c-e_1-e_2), (\alpha_1=e_1-e_2), (\alpha_j=e_j-e_{j+1}) for (2\le j<m), (\alpha_m=\sqrt2 e_m). The coefficients are (1,1,2,\ldots,2,\sqrt2).
* (\widetilde D_m): the preceding first roots and interior chain, with final fork roots (\alpha_{m-1}=e_{m-1}-e_m), (\alpha_m=e_{m-1}+e_m). The coefficients are (1,1,2,\ldots,2,1,1).

Deleting an interior node (i) leaves respectively

\[
B_i+B_{m-i},\qquad D_i+B_{m-i},\qquad D_i+D_{m-i}.
\]

After localization, all nonzero contacts are at interior nodes with common coefficient (d=1,2,2), respectively. The spatial coordinates of a normal with row (p) are, along the separating chain,

\[
z_k=-\sum_{\ell=k}^{m}p_\ell.
\]

For two rows supported within two nodes of (i<j), every coordinate indexed from (i+3) through (j-2) has normalized difference exactly (1/d) in absolute value. One row has no remaining mass there and the other has total contact mass (\tau/d). Therefore

\[
E\ge\frac{j-i-4}{d^2},\qquad
-(x_i,x_j)\ge\frac{j-i-4}{2}-2
\quad\text{if }j-i-4\ge4.
\]

Let (s) be the node of the narrow facet. Initially choose the deleted nodes as follows:

| Family | Initial (i,j) | Component rank requirements |
|---|---|---|
| (C_m) | (5,m-5) | both (B\) ranks at least 5 |
| (B_m) | (9,m-5) | (D\) rank at least 9, (B\) rank at least 5 |
| (D_m) | (9,m-9) | both (D\) ranks at least 9 |

If (s=i), increase (i) by one; if (s=j), decrease (j) by one. At the excluded ranks these alternatives are disjoint, so at most one move occurs. All retained-component thresholds remain true. The resulting negative pairings are at least

\[
\frac{m-19}{2},\qquad \frac{m-23}{2},\qquad \frac{m-27}{2},
\]

respectively. They reach 14 at (m=47,51,55).

## 4. The cycle energy for affine (A)

The affine (A_{n-1}) Gram matrix is the cycle Laplacian; all null coefficients equal 1. The normalized row difference has total charge zero, with positive mass 1 in the five-node interval about (i) and negative mass 1 in the interval about (j).

The Euclidean energy (E) is the electrical flow energy on this cycle. This follows either by solving the root-coordinate differences or from the Laplacian pseudoinverse. On the two long arcs separating the contact intervals, conservation makes the flow constant, with respective values (u) and (u-1). If their lengths are (L_1,L_2), all remaining edge energies are nonnegative, and

\[
E\ge L_1u^2+L_2(u-1)^2
\ge \frac{L_1L_2}{L_1+L_2}.
\]

The last inequality is the identity obtained by subtracting the right side and leaving (((L_1+L_2)u-L_2)^2/(L_1+L_2)). It covers arbitrary nonnegative contact distributions, not just single-node rows.

Choose opposite deleted nodes, avoiding (s) by rotating both one step if needed. Then

\[
L_1=\lfloor n/2\rfloor-4,\qquad L_2=\lceil n/2\rceil-4.
\]

For (n\ge136), both lengths are at least 64, so (E\ge32) and (-(x_i,x_j)\ge14). This recovers the linear separation mechanism of Nikulin's worked affine-(A) example with an explicit boundary.

## 5. Why these are genuinely framing walls of the narrow facet

This incidence bridge is essential. Avoiding (s) means that the narrow-facet root belongs to both retained edge supports. Its inverse diagonal is at least 1 in the path cases, and at least ((n-1)/n>1/2) in the cycle case. The Schur bound therefore gives

\[
|(x_i,\alpha_s)|<2.
\]

Thus their rank-two Gram is positive definite; the two hyperplanes meet **inside** hyperbolic space. More strongly, they share an actual ridge of (P). To see this, project a strict interior point (h) onto their common orthogonal complement. The inverse Gram is entrywise nonnegative, while the two pairings with (h) are negative. Consequently the projection has negative square and its pairing with every other wall remains strictly negative by acuteness. It is a point of the relative interior of that ridge. This argument does not rely on the terminating edge endpoint being finite.

Both terminating walls therefore belong to exactly the set controlled by the narrow-facet theorem. If their pairing is at most (-14), their hyperplane distance is at least (\operatorname{arccosh}7), contradicting the published strict upper bound. The four dimension bounds follow.

## Scope and the remaining problem

This proof closes the full-rank irreducible classical cusp calculation. The exceptional irreducible affine types have bounded rank at most 9; thus, using the standard affine classification, any **irreducible** cusp on the chosen narrow facet has ambient dimension at most 135. This is a conditional statement about that particular narrow facet. A noncompact polytope may have a narrow facet with only finite vertices, and a nonsimple ideal vertex has a product cusp diagram. Neither case is eliminated here.

For context, [Felikson–Tumarkin's primary theorem](https://arxiv.org/abs/math/0502413) already gives (n\le8) when the entire polytope is both simple and ideal. Their ideal hypothesis means every vertex is ideal, and is substantially stronger than the local hypothesis used here. The present bounds are proof-reconstruction lemmas, not a claim to improve that theorem or a verified new global bound for simple finite-volume polytopes.

For a product cusp the affine null vectors of the different factors are proportional but have independent relative scales. The row identity survives after choosing one common null vector. The localization and energy arguments above cannot simply be summed: an outgoing edge omits a node in **each** affine factor, and a terminating wall may distribute its contacts across those factors. Controlling that distribution and the narrow-facet incidence uniformly is the next unresolved problem. No claim that terminating walls contact only one component is used or justified here.

```sh
python scripts/fresh_nikulin_narrow_20261015.py --verify
```
