# Joint arithmetic witnesses and spectator-sensitive weights

The nine E8 faces in the sparse sealed-core weight obstruction are markings of **one compatible ten-wall local chamber**. Thus checking their joint arithmetic compatibility does not remove that obstruction. The D7 sparse witnesses behave differently: four require an isolated A1 in the full vertex diagram, provided the face remains unpaid by the existing H4 rule. This motivates weights depending on the full vertex forest.

The exact certificate is [the checker](../scripts/fresh_sqrt5_arithmetic_weight_obstructions_20261015.py) and [its data](fresh_sqrt5_arithmetic_weight_obstructions_20261015.json). Generation and saved replay pass. It uses the frozen independent Fraction arithmetic from the quadrilateral filter review, read-only; this is a new diagnostic, not an independent re-review of that arithmetic or of the inherited contact census.

## The coherent E8 chamber

Number the E8 tree by the chain 0–1–2–3–4–5–6 with node 7 attached to 2, and put \(\varphi=(1+\sqrt5)/2\). Add two norm-two roots with pairings \(-\varphi\) only at nodes 0 and 6 respectively. The two added roots are orthogonal. If \(J\) is the E8 inverse Gram matrix, their timelike residual squares have magnitudes

\[
d_0=\varphi^2J_{00}-2=4+2\sqrt5,
\qquad d_6=\varphi^2J_{66}-2=1+\sqrt5,
\qquad d_0=\varphi^2d_6.
\]

Since \(J_{06}=2\), their mutual product is \(\varphi^2J_{06}-\varphi d_6=0\). The certificate checks all 100 entries of each of nine literal marked-core identifications. The ten-root Gram matrix has inertia \((8,1,1)\); its nonidentity conjugate has inertia \((9,0,1)\). All entries are algebraic integers in \(\mathbb Q(\sqrt5)\), and the marked 5-edge square generates that field. These are the usual local arithmetic compatibility conditions, not a finite-covolume assertion.

An explicit nine-dimensional realization verifies a common strict timelike center and ten strict facet witnesses. All 45 eight-root supports are checked: the feasible cone has 21 finite and four spacelike extreme rays, with no null rays. A saved rational-in-the-field perturbation of a spacelike ray is strict on all ten walls and remains spacelike. Consequently this ten-wall chamber has infinite volume. No compact global extension or new high-dimensional example is claimed.

## Which sparse D7 witnesses require A1

Use the D7 chain 0–1–2–3–4–5 and the extra leaf 6 at 4. The complete six-witness sparse dual from the frozen component-only calculation has the following marked-pair orbits:

| Core classes | D7 marked pairs | Spectator condition |
|---|---|---|
| 535, 563 | \(\{1,4\},\{0,4\}\) | None |
| 843 | \(\{0,3\}\) | Isolated A1 |
| 967 | \(\{1,5\},\{1,6\}\) | Isolated A1 |
| 1032 | \(\{0,5\},\{0,6\}\) | Isolated A1 |
| 1052 | \(\{2,5\},\{2,6\}\) | Isolated A1 |

The checker reconstructs the literal corner forests, both D7 automorphisms, and all these pair orbits. For each of the six cores the frozen, independently reviewed contact census has no accepted nonempty **unpaid** external-white contact. A first external white node touching the core therefore either violates a necessary condition or creates an already-paid H4 configuration. If the face remains unpaid, every external white node is orthogonal to the core; each displayed isolated spectator is then an entire A1 component of the full vertex forest. The external forest may itself contain additional A1 components. This conditional sealing premise is inherited from the [contact report and review](fresh_sqrt5_quad_contacts_20261015.md).

The distance-two D7 baseline is 12. These six witnesses alone force total 14 without their A1-dependent demands and 21 with them. **The number 14 is not a valid all-face D7 bound.** Other sealed faces, in particular classes 785 and 860, still require corrections when no A1 exists. Their treatment belongs to the separate forest-sensitive rule; this package makes no assertion that the sparse six witnesses exhaust D7 faces.

The frozen component-only lower bound remains valid for weights that ignore the rest of the vertex diagram. Its sparse forcing faces cannot all be orthogonally inserted into the displayed D7³+E8 forest while preserving that exact forest: four of the D7 witnesses carry A1 spectators. This is a scope correction to that informal extension interpretation, not a change to the frozen linear-program certificate.

No dimension bound follows here. The exact joint E8 example defeats one proposed shortcut, while the D7 spectator condition identifies a concrete way around the component-only weight obstruction. Complete quadrilateral coverage remains a separate requirement.
