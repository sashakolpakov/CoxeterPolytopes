# Complete finite-endpoint domain of the new pi/9 edge

The actual open edge[58,59] in the [60-wall chamber](fresh_nonarithmetic_pi8_quad_completion_20261017.md) has positive support

\[
F=\{0,1,2,3,5,6,58,59\}=H_4+A_2+I_2(9).
\]

Its known finite endpoint also lies on55. There are exactly **seven possible finite-terminal contact profiles** for an additional wall: it is right to H4 and I2(9), and is either right to both A2 nodes or has a single bond3,4 or5 to node5 or6. Each profile has a complete one-parameter geometric endpoint domain. This certificate does not classify every finite-contact event against the remaining old walls, nor assert that several individually allowed normals are mutually compatible.

The [checker](../scripts/fresh_nonarithmetic_pi9_endpoint_domain_20261017.py) and [saved output](fresh_nonarithmetic_pi9_endpoint_domain_20261017.json) use exact recipes with outward rational intervals. Generation and fresh saved replay pass. A separate root-agent source/proof read accepts the profile reduction, conic coverage, tail monotonicity and neighboring-edge obstruction; it is not an independent arithmetic replay. All inputs are hash-bound, and no frozen construction is changed.

## Complete profile reduction

H4 has no nonright extension inside a finite diagram. The inverse diagonal of the I2(9) norm-two Gram is

\[
\frac{2}{4-4\cos^2(\pi/9)}>2.
\]

Any nonzero contact has magnitude at least1 and therefore already exceeds the available norm-two Schur budget. Both I2(9) contacts must be right.

The A2 inverse is one third of the matrix with diagonal2 and off-diagonal1. A single bond of magnitude c costs2c²/3. Labels3,4,5 give cost below2; label6 gives equality and larger labels give greater cost. Two nonzero contacts have cost at least2. Strict positive definiteness of a finite endpoint leaves exactly the seven profiles above. This excludes all larger labels analytically; it is not a label cutoff.

## One surviving hyperbola branch

Normalize the known finite endpoint to a future unit timelike vector k0. The old wall55 is orthogonal to F and to this endpoint, so k1=n55/sqrt2 is a unit spacelike vector in the same Lorentz2 complement. Put

\[
A=(k_0-k_1)/\sqrt2,\qquad B=(-k_0-k_1)/\sqrt2.
\]

For each profile the exact positive-span projection p and residual square rho=2−p² are saved. All possible normals in the surviving orientation are

\[
r(t)=p+\sqrt{\rho/2}\,(tA+B/t),\qquad t>1.
\]

The other orientation has positive product with55 for every t>0, so fails acuteness. In the surviving branch the corresponding future endpoint is

\[
x(t)=tA-B/t,\qquad x(t)^2=-2,\qquad x(t)\cdot r(t)=0.
\]

Its55 product is−t+1/t, forcing t≥1. At t=1 it is the existing finite vertex. A distinct acute added wall through that vertex cannot share a strict interior: it lies in the positive-definite linear span of the nine incident normals, and inverse positivity makes its coefficients nonpositive, hence its product with every old strict interior point positive. Thus t=1 is excluded.

All other old endpoint inequalities hold strictly for **every** t>1. The output verifies A·n_i<0 outside F and (A−B)·n_i<0 outside F∪{55}. Therefore

\[
x(t)\cdot n_i=(A\cdot n_i)(t-1/t)+((A-B)\cdot n_i)/t<0.
\]

Every r(t) with t>1 is strict at the old finite endpoint. Whenever its native old contacts are acute, a sufficiently small perturbation of that endpoint toward the old strict center gives a common future interior. The same argument applies to a finite mutually acute family; mutual Coxeter compatibility must still be checked separately.

## Certified coherent tails

Every native old contact is saved in the form

\[
p\cdot n_i+\sqrt{\rho/2}\,(a_i t+b_i/t).
\]

At the following dyadic thresholds all non-support contacts are strictly below−2, with strictly negative derivatives:

| A2 profile | Valid all-dotted tail |
|---|---|
| zero | t≥4 |
| one bond3 at5 or6 | t≥4 |
| one bond4 at5 or6 | t≥4 |
| one bond5 at5 or6 | t≥8 |

These conclusions hold throughout the indicated real intervals. Since a_i<0, the derivative a_i−b_i/t² is negative for every larger t if b_i≥0; if b_i<0 it decreases toward a_i, so its threshold sign also persists. Thus these are actual continuous families of individually coherent finite-endpoint walls, not numerical samples. Their simultaneous use is not asserted.

## A single wall cannot also terminate[30,59]

There is a useful complete obstruction to collapsing the next coordination problem into one wall. The support of the neighboring edge[30,59] is H4+A2+I2(5). A proposed common terminal is already right to59. A nonright contact to30 costs at least

\[
d_5=\frac{2}{4-\varphi^2}.
\]

Any nonzero A2 contact adds at least2/3, and d5+2/3>2. Hence a nonright30 contact forces the zero A2 profile. For that profile only labels2 and3 to30 fit the I2(5) budget: label4 already costs2d5>2. Label2 is t=1 and is excluded above. Label3 gives a unique positive quadratic root

\[
t\approx1.06862520816.
\]

At that point the old49 product is strictly positive, approximately1.86641852964. It therefore fails acuteness. The exact quadratic and positive interval are saved. **No single acute wall can provide finite terminal vertices for both[58,59] and[30,59], at any field or finite labels.** Two coordinated walls remain possible.

This package supplies finite-terminal domains and an all-label single-wall obstruction. It establishes neither compactness of the whole target H3 face nor finite volume of the ambient chamber. The strict tails demonstrate that the remaining edge itself is viable; they do not eliminate its neighboring closure obligations.

The finite-terminal requirement is also necessary in every finite-volume completion retaining this pi/9 edge. By the standard classification of irreducible Euclidean Coxeter diagrams, their finite nonright labels are3,4,6. A finite-volume ideal vertex link is a product of these affine components. Its incident roots therefore cannot contain a pi/9 pair. Thus the other endpoint of the retained edge cannot be ideal; no simplicity assumption at the cusp is needed.
