# Compact dimension 21 is reduced to two fields

Combining the [new eight-field theorem](fresh_arithmetic_high_label_20261016.md)
and its [independent review](fresh_arithmetic_high_label_20261016_review.md)
with the previous Q(sqrt2) and K7 proofs leaves precisely

\[
\boxed{\mathbb Q(\sqrt5),\qquad\mathbb Q(\sqrt2,\sqrt5)}
\]

unresolved for a compact arithmetic Coxeter polytope in dimension 21.
This is a logical composition of internally reviewed proof candidates,
not an exclusion of dimension 21. The all-arithmetic upper-bound candidate
remains 21; no compact example in that dimension is asserted.

The [binder](../scripts/fresh_compact21_field_reduction_20261016.py) and
[certificate](fresh_compact21_field_reduction_20261016.json) account for
all thirteen branches in the published dimension-at-least-14 field list.
The rational branch is excluded separately: a rational indefinite form of
rank at least five is isotropic by Meyer, precluding cocompactness of the
corresponding simplest-type arithmetic lattice. The twelve nonrational
branches consist of eight new bound-20 candidates, the two prior bound-20
candidates, and the two explicitly open fields above.

For Q(sqrt2), the binder includes both the all-size quadrilateral closure
review and its separately reviewed weight calculation. The latter's
conditional coverage premise is discharged by the former. Its dimension-21
entry has strict lower weight 40 and vertex upper weight 40. The cubic-seven
review supplies its complete geometric and counting proof. The new eight
fields have independent cyclotomic automorphism certificates, the complete
local support reduction and the same strict dimension threshold.

The field list is from Vinberg's Theorem 3, also restated in
[Nikulin's introduction and Section 3](https://arxiv.org/pdf/0708.3991).
The three-field restriction that starts at dimension 22 is not used to
remove a dimension-21 branch. A separate read-only scope review checked
this distinction, the rational branch, and the Q(sqrt2) dependency pairing.

Generation and saved replay check 96 recursive manifests. These hash
checks preserve the exact inherited inputs; they do not rerun their large
enumerations or replace mathematical validation.

```sh
python3 scripts/fresh_compact21_field_reduction_20261016.py
python3 scripts/fresh_compact21_field_reduction_20261016.py --verify
```

Neither remaining field is covered by the new general field hypothesis.
The F=Q(sqrt2,sqrt5) branch has only order-two/three arithmetic Lanner
diagrams by the published arithmetic Lanner classification, but triangles
such as (2,5,8) and (4,5,5) defeat the present distinguished-bond argument.
The complete Q(sqrt5) weighting currently excludes dimensions 22 upward,
not 21. Both branches require new work.
