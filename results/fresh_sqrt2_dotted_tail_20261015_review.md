# Independent review of the C5 pendant-tail lemma

**PASS.** The [independent checker](../scripts/fresh_sqrt2_dotted_tail_20261015_review.py) and [certificate](fresh_sqrt2_dotted_tail_20261015_review.json) verify the [primary classification](fresh_sqrt2_dotted_tail_20261015.md), including its all-attachment argument within the stated topology. No gap was found. This is a local classification, not a complete quadrilateral theorem or a dimension bound.

The checker imports neither the producer nor a symbolic library. It uses a separate `Fraction` implementation of \(\mathbb Q(\sqrt2)\), permutation cycle-cover determinant expansion, congruence pivots, and exact residual geometry. The source/data bindings and every primary manifest entry are checked.

## The all-attachment argument

Use the primary names: opposite pairs \(a,y\) and \(b,d\), white chain \(u-v-w\), the C5 \(b-u-v-w-d-b\) with its sole 4-bond \(b-d\), and contact \(a-v\) of label 3. In the double family \(y-v\) also has label 3. All other unlisted core pairings are right. Every extra connected white tree is attached by one edge to one core node; an orthogonal white component is harmless.

Both corners containing \(a\) have a D5 principal diagram. Its connected spherical extensions are simply laced D/E diagrams with only one branching vertex. Consequently a pendant contact at \(v\) creates a forbidden fourth arm, while a contact at \(u\) or \(w\) creates a second branching vertex in one of the two corners. At \(a,b,d\), the only remaining possibility is one all-3 path. In particular, the first contact at \(b\) or \(d\) is **exactly** label 3, not merely at least 3. Together with the opposite 4-bond and the existing white contact it then gives the literal affine fork with arm labels 4,3,3. Compactness excludes it. The independent checker verifies its positive null vector and inertia.

The path at \(a\) gives D5,E6,E7,E8 for lengths 0,1,2,3; length 4 gives affine E8. This is checked by independently reconstructed matrices. Longer paths contain the forbidden prefix, and branching or a second attached path is already excluded by the finite-component classification. The double family gives the same restriction at \(y\).

In the single family the \(y\)-component is separate from the C5 part in each finite corner. Its reduced unit diagonal \(q\) is positive and at most one. It is rational unless it is the isolated I2(8) component, which gives \(q=(2-\sqrt2)/4\). Indeed every other permitted finite component has rational squared bonds, so exact leaf elimination remains rational; a label at least 6 cannot have a nonright neighboring edge in a finite diagram. Multiple trees at \(y\) are covered by treating their union with \(y\) as its one finite component. A nonzero attachment has strictly positive projection cost, so \(q=1\) forces no attachment.

These arguments cover arbitrarily large permitted pendant trees without enumerating their size. They do not cover an external component with two core attachments, extra core edges, shared core vertices, or changed core labels.

## Independent rank and integrality checks

Permutation expansion of the seven-by-seven determinant reproduces both identities, with \(p,q\) the reduced unit diagonals and original norm-two contact \(-r\):

\[
\det G_{\rm single}=2(1+\sqrt2)\{r^2-q(4p-1+\sqrt2)\},
\]

\[
\det G_{\rm double}=2(1+\sqrt2)
\{r^2+(1-\sqrt2)r-4pq+(1-\sqrt2)(p+q)\}.
\]

Independent inversion of each finite path gives \(p=(\ell+2)/(2(\ell+1))\). Corank one is imposed on the **complete** affected star after all allowed extra white components have been eliminated. Orthogonal components contribute positive factors; a determinant calculation on an arbitrary smaller subdiagram would not suffice.

In the single family, the original two-cycle is

\[
D=r^2=q\left(\frac{\ell+3}{\ell+1}+\sqrt2\right).
\]

If \(q\) is rational, integrality in \(\mathbb Z[\sqrt2]\) makes its \(\sqrt2\)-coefficient \(q\) an integer. Thus \(q=1\), and the rational coefficient permits only \(\ell=0,1\). All four exceptional I2(8) values are independently found nonintegral. The survivors have \(D=3+\sqrt2\) and \(2+\sqrt2\), respectively.

For the double family, the original triangle \(a-v-y\) has two contacts of magnitude one. Its cyclic product is \(-r\), so integrality requires \(r\in\mathbb Z[\sqrt2]\), a stronger requirement than integrality of \(r^2\). Writing \(r=A+B\sqrt2\) with integers \(A,B\), the rank equation's \(\sqrt2\)-coefficient gives

\[
2AB+B-A=p+q.
\]

Among \(p,q\in\{1,3/4,2/3,5/8\}\), only \(p=q=1\) has integer sum. The checker evaluates all 16 pairs. At that pair the equation factors as \((r-1-\sqrt2)(r+2)\), leaving the unique positive root \(r=1+\sqrt2\). This independently proves the integrality classification without using the producer's field-square-root routine. The review also checks every non-null contact recorded by that routine satisfies the equation and its recorded integrality status; it does not separately reproduce the intermediate nonintegral square-root membership decisions, which are unnecessary to this alternative proof.

The distinction between original and computational normalization matters. For single-family signature checks the producer replaces the original \(y\) by \(\sqrt D\,y\), giving square \(2D\) and contact \(-D\). Since \(D,D'>0\), this is a valid congruence at both relevant embeddings. The arithmetic exclusions above concern the original norm-two two-cycle \(D\), not cyclic products formed after that computational rescaling. The double-family cycle condition similarly uses the original unit-magnitude contacts.

## Signatures, compact local faces, and remaining scope

For all three surviving matrices the independent checker reconstructs the physical and conjugate inertia, all four finite corners, and all **509** nonempty principal-subset tests excluding positive-semidefinite singular blocks. It additionally eliminates the white support, solves for all four residual vertex rays, and verifies negative norms, negative pairings with both nonincident sides, and negative pairings between all vertex rays. Thus the four finite corners are feasible and future-oriented in a common cone; their compact quadrilateral is not inferred from rank alone.

The indicated D5 marked long-leaf/short-leaf pairs and E6 marked near/far pairs agree with these literal corner graphs. Assigning half weights to the stated two pairs of an **entire** D5 or E6 component pays the three local survivors as described in the primary. This does not establish the corresponding inequality for other quadrilateral cores. External white components with multiple core attachments can produce Schur off-diagonal corrections and remain outside the present lemma.

Generation and saved replay pass:

```sh
python3 scripts/fresh_sqrt2_dotted_tail_20261015_review.py
python3 scripts/fresh_sqrt2_dotted_tail_20261015_review.py --verify
```

The complete result is exactly the primary's two fixed C5 attachment families with permitted pendant trees and orthogonal components. No global dimension improvement, complete ambient polytope, exhaustive quadrilateral classification, or novelty claim is certified.
