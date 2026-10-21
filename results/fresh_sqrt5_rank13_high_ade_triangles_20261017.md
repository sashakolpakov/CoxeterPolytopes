# Triangles through the joint high-baseline ADE set

The [joint rank13 endpoint census](fresh_sqrt5_rank13_high_ade_endpoints_20261017.md) contains215 ADE forests with baseline at least12. This package completes the necessary local arithmetic and finite third-corner tests for triangles with two vertices in that set, allowing the two forests to differ.

There are24,523 accepted literal occurrences:21,820 with framing label5 and2,703 with framing label10. The checker tests232,772 endpoint/label choices and1,084,212 third corners. The possible root-pair counts range from9 to72 per selected vertex type. D13 has9, A13 has10, and D5+2D4 has19 once transitions to other selected forests are allowed. The maximum72 occurs at A1+4A3, whose baseline is12.

These counts are upper lists of individually admissible local corners. They do not assert simultaneous incidence or compact ambient realization. They show that a uniform half-unit payment at every possible two-special-triangle corner is much too coarse for some joint-set types; for A1+4A3 that estimate is36 versus the target12. This is a limitation of that estimate, not a proof that a refined global incidence argument fails.

## Exact arithmetic reduction

Let J be the inverse of the old positive-definite ADE Gram. If i is the removed root and S the new root's retained contacts, write

`a=J_ii`, `b=−sum_(j in S) J_ij`, `c=2−sum_(j,k in S) J_jk`.

With framing product g, the new Schur entry is `c−2bg−ag²`. The endpoint graph has only zero and minus-one retained products. The complete finite-label field restriction leaves2,3,4,5,6,10.

Labels2 and3 make the whole edge Gram rational, so its physical negative direction contradicts conjugate positivity. For labels4 and6, a fixed label-three path between the framing endpoints violates the cycle-field condition. Without such a path, b=0 and g² is rational, making the Schur signs identical at the two embeddings. These labels are therefore excluded as well.

For label5, g=−phi and the checker evaluates the Schur entry and its conjugate with exact Fraction arithmetic in Q(sqrt5). For label10, the framing edge must be a bridge, b=0 and g²=2+phi; the same exact sign test applies. Physical negativity and strict conjugate positivity are required. The remaining old Gram is positive definite at both embeddings, so these tests give the full edge inertia.

After the other released root is deleted, only the component containing the framing bond can fail to be an induced finite ADE subdiagram. A label-five component must be a path with that bond at one end and at most four roots: I2(5), H3 or H4. A label-ten component must be I2(10). This exact finite classification is applied to every other released root. The output retains the endpoint attachment index, destination type, framing label, released root and exceptional third component for every accepted case.

## Why three selected vertices cannot form a triangle

Choose K-rational normals with totally positive norms a_i and adjoin the normalizing square roots `sqrt(2/a_i)`. Every nonphysical embedding of K extends to a real embedding of this extension, because each conjugate a_i is positive. If all three triangle corners are ADE, the full normalized star Gram is rational and is fixed entrywise. Its physical negative direction would then occur in the positive-definite conjugate form, a contradiction. The forests at the three vertices need not agree.

The [source](../scripts/fresh_sqrt5_rank13_high_ade_triangles_20261017.py) and [saved data](fresh_sqrt5_rank13_high_ade_triangles_20261017.json) pass generation and a fresh saved replay:

```sh
python scripts/fresh_sqrt5_rank13_high_ade_triangles_20261017.py
python scripts/fresh_sqrt5_rank13_high_ade_triangles_20261017.py --verify
```

The graph inputs are inherited from the complete endpoint census; the arithmetic is a direct Schur calculation using the frozen standalone Fraction field engine. Joint all-special quadrilateral classification, more precise payments at the third vertices, and simultaneous face incidence remain unresolved. No joint dimension bound follows from this package.
