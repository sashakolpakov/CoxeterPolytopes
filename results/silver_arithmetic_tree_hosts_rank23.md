# Silver arithmetic tree hosts in rank 23

This exact census separates two easily confused notions.  The product
\(c=2+\sqrt2<4\) is the finite Coxeter label-eight product
\(4\cos^2(\pi/8)\), although it generates \(\mathbf Q(\sqrt2)\).  A genuine
dotted edge instead has product greater than four; the explicit mild choice
below is \(p=3+2\sqrt2\).

## One field-generating label-eight leaf

A leaf with Cartan pair \((-c,-1)\) was attached at every labelled node of
the finite skeletons \(A_{22},D_{22},B_{22}\).  The resulting symmetric Gram
is integral over \(\mathbf Z[\sqrt2]\).  Exact two-place inertia and all 23
inverse-Gram diagonal entries give:

| skeleton | attachments with inertia `(22,1)` and conjugate PD | best attachment | finite / ideal / hyperideal |
|---|---|---:|---:|
| `A22` | `0,1,2,3,18,19,20,21` | endpoint `0` or `21` | **3 / 0 / 20** |
| `D22` | long-arm nodes `19,20,21` | terminal `21` | **3 / 0 / 20** |
| `B22` | long-arm nodes `19,20,21` | terminal `21` | **3 / 0 / 20** |

Thus attaching at the branching or minuscule end does not beat the terminal
\(A_{22}\) seed.  The complete labelled scan contains 66 attachments; cases
not listed fail conjugate positive definiteness (selected inertia is still
Lorentzian).

For the terminal \(D_{22}\) host, the inverse diagonal begins
\(19-\sqrt2,(21-\sqrt2)/4,(21-\sqrt2)/4,18-\sqrt2,\ldots\) and ends
\(1-\sqrt2,-\sqrt2,(1-\sqrt2)/2\).  Exactly the last two skeleton rays and
the new leaf ray are finite.  The terminal \(B_{22}\) host has the same sign
pattern; its last three entries are
\((1-\sqrt2)/2,-\sqrt2/2,(1-\sqrt2)/4\).

## Structural obstruction for two genuine dotted edges

Deleting one wall gives a simplex vertex link.  If its principal Gram still
contains a genuine dotted edge, the corresponding 2-by-2 principal block has
negative determinant, so that ray is hyperideal.  Hence a deletion can be
finite or ideal only when its wall is incident to **every** dotted edge.
Two distinct adjacent dotted edges have only one common endpoint and force at
least \(R-1\) hyperideal rays in rank \(R\); two disjoint dotted edges force
all \(R\).  In ranks 23--30 the adjacent lower bounds are therefore 22--29.

The bound is sharp within arithmetic silver trees.  In rank 23, put dotted
products \(p=3+2\sqrt2\) on `0-1` and `0-2`, and label-three edges on the arm
`0-3-...-22`.  With symmetrizer
`(1,p,p,1,...,1)`, the selected Schur complement is

\[
  2-p-20/21=22/21-p=-41/21-2\sqrt2<0,
\]

while its conjugate is positive.  Thus the exact inertias are `(22,1)` and
`(23,0)`, the minimal field is \(\mathbf Q(\sqrt2)\), and the cone has exactly
**1 finite, 0 ideal, 22 hyperideal** rays.

This is a cone-triage theorem, not a compact example or a Vinberg
nonreflectivity theorem.  It shows that genuine multi-dotted simplicial seeds
cannot supply the hoped-for small-hyperideal starting cone; the best tested
finite-skeleton seed remains the one-label-eight terminal attachment.
