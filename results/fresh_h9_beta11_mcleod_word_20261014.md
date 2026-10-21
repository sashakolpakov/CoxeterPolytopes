# The transported McLeod symmetry is already a reflection word

**PASS.** The integral infinite-order matrix proposed in the frozen
[literature comparison](fresh_h9_beta11_literature_20261014.md) belongs to
the known reflection subgroup of the determinant-11 lattice. Consequently
this particular matrix supplies no nontrivial chamber-stabilizer
obstruction. This result does not decide reflectivity of the lattice.

The [standalone checker](../scripts/fresh_h9_beta11_mcleod_word_20261014.py)
and [certificate](fresh_h9_beta11_mcleod_word_20261014.json) use fresh
Fraction-pair arithmetic in \(\mathbb Q(\phi)\), with no producer or
project arithmetic imports. They bind the frozen literal McLeod transport
and the [height-100 both-norm prefix](fresh_h9_beta11_vinberg_20261014.json).
Coordinates are \((x,E8_0,\ldots,E8_7,t)\), and coefficient pairs mean
\(a+b\phi\).

For a saved root \(r_i\), the checker constructs

\[
 S_i=I-\frac{2}{r_i^TGr_i}\,r_i(r_i^TG).
\]

Every reflection used below has integral entries, squares to the identity,
preserves G, and negates its positive-norm root. The eleven used indices
are \(0,1,\ldots,9,11\). Root 9 has norm \(2\beta\); the others have
norm two. Thus the new norm class is present in the word.

The literal chronological list of left multiplications is

```text
2,0,3,2,4,3,1,5,4,3,2,0,6,5,4,3,1,2,3,4,5,7,6,5,
4,3,1,2,0,3,2,4,3,1,5,4,3,2,0,8,9,8,9,8,9,8,9,11,
8,9,8,11,8,9,0,2,3,1,4,3,2,0,5,4,3,1,2,3,4,5,6,5,
4,3,1,2,0,3,2,4,3,1,5,4,3,2,0,7,6,5,4,3,2,8
```

Writing its entries as \(i_1,\ldots,i_{94}\), the exact identity is

\[
 S_{i_{94}}\cdots S_{i_1}g=I,
 \qquad\text{equivalently}\qquad
 g=S_{i_1}\cdots S_{i_{94}}.
\]

Both multiplication directions are checked independently. After the first
54 left multiplications, the residual fixes the Vinberg basepoint
\((2,0^8,1)\), but its matrix is not yet the identity. After 94, the whole
matrix is exactly the identity. The intermediate matrix is serialized.
The transport remains infinite order: its real trace is greater than ten,
which excludes finite projective order. Infinite order alone therefore
does not distinguish a chamber symmetry from an element of the reflection
group.

The reflection group acts freely on its chambers, so this nonidentity
element cannot stabilize a chamber of its full mirror arrangement. Later
roots cannot change the established subgroup membership. No new root
enumeration or completeness claim is needed for the word certificate;
the known integral reflections already suffice.

Generation, saved replay, syntax, source manifests and report links pass.
No frozen predecessor was modified.

```sh
python scripts/fresh_h9_beta11_mcleod_word_20261014.py --verify
```
