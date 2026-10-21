# Review of the arbitrary-normal H20 finite-angle restriction

The [main proof](fresh_h20_crystallographic_angle_gate_20261014.md) was reviewed and its primary saved-result replay passed. No mathematical gap was found in the stated conclusion: an actual finite ridge of a finite-volume Coxeter H20 hyperplane section of the literal 210-wall source can have only labels 2, 3, 4, or 6. This is an angle restriction, not an arithmeticity or general nonexistence theorem.

The additional [separate Fraction checker](../scripts/fresh_h20_crystallographic_angle_gate_20261014_review.py) and [certificate](fresh_h20_crystallographic_angle_gate_20261014_review.json) reconstruct the six new interval tables without importing the primary arithmetic. In each table, the squared-pairing excess over the norm product is actually linear in `q`, since its quadratic coefficient cancels. Rational endpoint signs therefore certify the whole open interval. The check also reproduces the label-5 long-A2 inequalities and the rational/algebraic endpoints used to exclude all B2 finite labels. It does not repeat the previous H3 irreducibility enumeration or the entire old spherical-triple catalog in a second implementation; those remain explicit dependencies of the reviewed proof and primary replay.

The geometry and scope checks were as follows.

- A projected positive-definite triple has positive-definite old Gram, since the old Gram adds a positive semidefinite rank-one term. Thus the literal no-triangle catalog applies to every facet of an assumed actual dihedral face, with no restriction on the section normal's other contacts.
- For `m >= 5`, `m != 6`, an ideal point on the face would put the label-m pair into a bounded Euclidean Coxeter link, which has no such affine diagram. The face is compact. For `m >= 7` there is no connected finite extension; for `m = 5` the previously proved all-H3 obstruction removes both possible extension orientations. All remaining boundaries are orthogonal to the dihedral pair.
- Two nonzero contacts on an old nonorthogonal pair force any modified boundary to create a forbidden old positive-definite triangle. For an old orthogonal pair, the literal catalog leaves exactly the three contact ratios used in the proof. The modified-wall pair tables then leave only unchanged perpendicular boundary normals on the required smaller compact face.
- With one nonzero contact, the four old pair types are exhaustive. The preceding short-A2 proof transfers to an assumed actual dihedral ridge: its original source-edge hypothesis was used to establish that ridge's actuality, not in the subsequent boundary equations or rational compact-face argument. The long-A2 exceptional `(m,n)=(8,8)` step and the two B2 cases require no additional vanished contacts.
- A compact face bounded solely by unchanged source roots has those roots spanning its full Lorentzian space. This makes the space rational, and the old integral reflection matrices preserve its full intersection lattice. The resulting uniform reflection lattice would have finite index in the integral orthogonal group, contradicting rational isotropy in the stated dimensions. Rational coordinates alone would not justify this step; the spanning and lattice-preservation premises are both present.

The proof assumes the original source inequalities and an actual finite-volume Coxeter section. It does not apply automatically after deleting source facets, subdividing a non-Coxeter section, or changing the ambient chamber. It also places no restriction on dotted distances or cyclic-product integrality.

```sh
python scripts/fresh_h20_crystallographic_angle_gate_20261014.py --verify
python scripts/fresh_h20_crystallographic_angle_gate_20261014_review.py --verify
```
