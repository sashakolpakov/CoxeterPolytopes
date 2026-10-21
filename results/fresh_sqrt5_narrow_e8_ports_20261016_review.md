# Independent review of the central-E8 narrow-facet restriction

**PASS:** no minimality-14 edge at a central E8 deletion can return from
`D5+2E8` to that same vertex type. The complete finite endpoint census has
29 types, and each of the 28 other types has at most four matching edge
ports. Consequently `D5+2E8` vertices occupy at most four fifths of the
vertices of a selected narrow facet in compact arithmetic dimension 21
over `Q(sqrt5)`. This does not exclude that dimension.

The [reviewer](../scripts/fresh_sqrt5_narrow_e8_ports_20261016_review.py) and
[certificate](fresh_sqrt5_narrow_e8_ports_20261016_review.json) use independent
quadratic arithmetic and the previously reviewed finite-diagram generator.
They do not import the [primary producer](../scripts/fresh_sqrt5_narrow_e8_ports_20261016.py).
The review covers the central-port argument and facet incidence conclusion;
the producer's other copied-port sharp bounds and example edges are outside
this review.

Deleting the central E8 root leaves arms `A1+A2+A4`. To restore E8 in the
common support `A1+A2+A4+D5+E8`, the new central root must attach to the A1
and one endpoint of each longer arm. These are exactly four possibilities;
the intact D5 and E8 components cannot be enlarged into another required
component. Independent inverses of the three rational chain Grams give
cross costs `59/30,49/30,41/30,31/30`, with residual norm `1/30` for each
framing root.

The shared A1 has product minus one with both framing roots. Their
triangular cyclic product therefore equals the unknown framing product
`g`, so arithmetic cyclic integrality gives `g in Z[phi]`. Positive
definiteness at the other embedding requires `q-1/30<g'<q+1/30`.
Minimality 14 gives `-14<g<=0`. From `g-g'=b*sqrt5` one gets
`-7<=b<=-1`, and the safely enlarged interval `-5<=a<=1` covers all
possibilities for `g=a+b*phi`. Direct exact enumeration leaves only
`g=1-phi` at cost `49/30`; it lies strictly between minus one and zero,
so is neither a finite nonright Coxeter contact nor a parallel/dotted one.

Independently deleting every root from each of the 68 finite connected
diagram types of rank at most 21 and appending untouched common-support
components gives all 29 endpoint types. A second deletion count gives their
full matching-port multiplicities. It reproduces the entire primary table
after identifying the usual `A2/I2(3)` and `B2/I2(4)` aliases. The maximum
outside `D5+2E8` is four, attained precisely by
`2A2+A4+D5+E8`, `2A4+D5+E8`, and `A1+A2+2D5+E8`.

On the selected narrow facet, let `Dc` count source vertices whose facet
normal is a central E8 root, `Do` the remaining source vertices, and `O`
all other vertices. At least one central deletion edge lies in the facet
at every source vertex, and two do unless the facet normal is central.
Every such edge has minimality 14 and hence ends at another vertex type.
Counting these edges at their other endpoints gives

`Dc+2Do <= sum_U m(U)V_U <= 4O`.

Thus `D= Dc+Do <=4O`, or `D/(D+O)<=4/5`. The existence of the selected
narrow facet, and minimality 14 for its edge chambers, are the published
inputs in [Nikulin, Section 4, printed pages 10–12](https://arxiv.org/pdf/0708.3991).
No claim is made about whole-polytope frequencies or about a chosen small
two-face containing one of the central deletion edges.

Generation and fresh saved replay passed. This is a genuine additional
incidence restriction, but it needs a further global argument to produce
a lower dimension ceiling.

```sh
python3 scripts/fresh_sqrt5_narrow_e8_ports_20261016_review.py
python3 scripts/fresh_sqrt5_narrow_e8_ports_20261016_review.py --verify
```
