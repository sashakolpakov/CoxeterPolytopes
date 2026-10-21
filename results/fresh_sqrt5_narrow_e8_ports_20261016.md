# A central E8 return obstruction on a narrow facet

**New scoped incidence restriction:** in a compact arithmetic Coxeter polytope of dimension 21 over Q(√5), vertices of type **D5 + 2E8** occupy at most **4/5 of the vertices of a selected minimality-14 facet**. More precisely, the weighted edge-count inequality below distinguishes the facet normal's position in the vertex diagram.

This does not exclude dimension 21 or E8 vertex components. The [exact checker](../scripts/fresh_sqrt5_narrow_e8_ports_20261016.py) and [certificate](fresh_sqrt5_narrow_e8_ports_20261016.json) first test copied attachments, then exhaust all same-type returns at the three relevant E8 ports. These are different assertions: some noncopied returns survive even when copied returns fail.

## Exact scalar test

Use norm-two roots, `phi²=phi+1`, and the sign convention that distinct wall products are nonpositive. For a copied root of a finite E8 Gram `A`, put

\[
c=(A^{-1})_{ii},\qquad s=1/c,\qquad q=2-s.
\]

If the framing product is `g`, eliminating the common roots gives `[[s,g−q],[g−q,s]]`. Conjugate positive definiteness requires

\[
q-s<g'<q+s.
\]

The eight inverse diagonals are `4,14,30,20,12,6,2,8`. For the three values **14, 20 and 30**, no copied-root arithmetic edge has `−14<g≤0`. The sharp smallest possible absolute product is

\[
|g|=3+8\phi>14.
\]

Actual finite local edges attaining this value are checked, including their signatures, positive conjugates and finite endpoint rays. Thus this is a sharp obstruction for copied attachments, not a numerical cutoff.

The finite search has a proved box. Whenever `g=a+b phi` is integral, `−14<g≤0` and `0<g'<2`,

\[
b=(g-g')/\sqrt5\quad\Longrightarrow\quad -7\le b\le-1,
\qquad -4\le a\le1.
\]

All 42 pairs are checked exactly. The sharp lower-bound check uses the equally justified enlarged range `−8≤b≤−1` below `|g|=3+8phi`.

Integrality is not inferred from normalization alone. For the central reattachments below, the common A1 root has product −1 with both framing roots, so the corresponding triangular cyclic product is `g`. The arithmetic cyclic-product criterion therefore puts `g` in the full integer ring `Z[phi]`. The nondegenerate edge Gram is positive definite at the other embedding. These are the relevant arithmetic inputs from [Nikulin, §2](https://arxiv.org/pdf/0708.3991).

## All central same-type returns fail minimality 14

Deleting the central E8 root, where `c=30`, leaves **A1 + A2 + A4** in that component. In D5 + 2E8, the complete common edge diagram is therefore

\[
S=A1+A2+A4+D5+E8.
\]

To return to D5 + 2E8, the new root must restore E8 from the first three components. The checker exhausts all subsets of the 20 common roots of sizes zero through three, using exact tree isomorphism to retain the desired endpoint type. This covers every possible endpoint: its new bonds must be 3 and its degree at most three. Exactly **four** central reattachments survive. They select either endpoint of A2 and either endpoint of A4, and all retain the common A1 contact.

Their cross-projection costs and conjugate intervals are:

| Cross cost `q` | Required interval for `g'` | Narrow integral candidates |
|---|---|---|
| 59/30 | (29/15, 2) | None |
| 49/30 | (8/5, 5/3) | `g=1−phi`, rejected |
| 41/30 | (4/3, 7/5) | None |
| 31/30 | (1, 16/15) | None |

The one scalar candidate has `−1<g<0`. It cannot be a Coxeter product: a nonright finite angle has norm-two product at most −1, while a dotted pair has product below −2. Hence **no minimality-14 central-E8 edge can return to D5 + 2E8**, even with a different attachment to the common diagram.

The distinction from copied attachments matters. The other two obstructed copied ports have valid narrow noncopied returns. At `c=14`, one uses `g=−2−5phi`. At `c=20`, returns use `g=−1−4phi` or `g=−2−5phi`. The checker verifies actual finite endpoints and conjugate positive definiteness for all three examples. No central-port conclusion is silently applied to those ports.

## From an edge obstruction to a facet incidence inequality

Nikulin's narrow-facet construction supplies a facet `F` whose framing hyperplanes satisfy minimality 14. An edge lying in `F` has an edge chamber of minimality 14. His simple-polytope estimate also supplies a triangular or quadrilateral two-face of `F` in the dimensions relevant here. It does not require that this small face meet a specified E8 port. [Nikulin, §4, especially printed pp.10–12 and §4.2](https://arxiv.org/pdf/0708.3991).

At a D5 + 2E8 vertex on `F`, there are two central E8 root ports. The facet normal is at most one of them. Deleting any other central port gives an edge contained in `F`; its opposite endpoint must have a different vertex type by the preceding result. Let:

- `D_c` count these exceptional vertices where the facet normal is a central E8 root;
- `D_o` count the other D5 + 2E8 vertices on `F`;
- `V_U` count other vertices on `F` of type `U`.

The number of distinguished outgoing edges is `D_c+2D_o`. Each has common diagram `S` above. The complete finite-diagram census gives **29** one-root completion types of `S`. Removing the forbidden D5 + 2E8 return leaves 28 types. For each type, the checker counts every root deletion yielding `S`; call this number `m(U)`. Every other vertex receives at most `m(U)` distinguished edges, and **m(U)≤4**. The maximum occurs precisely for

\[
2A2+A4+D5+E8,\qquad 2A4+D5+E8,\qquad A1+A2+2D5+E8.
\]

Counting each distinguished edge at its nonexceptional endpoint gives

\[
D_c+2D_o\le\sum_U m(U)V_U\le4O,
\qquad O=\sum_U V_U.
\]

Consequently, with `D=D_c+D_o`,

\[
\frac{D}{D+O}\le\frac45.
\]

The finite endpoint census uses the complete 68-type rank-21 catalogue already used by the frozen arithmetic models. Its completion argument is exhaustive: the new root lies in one connected finite component; deleting it gives a submultiset of `S`, and every other component stays unchanged. The certificate saves all 29 types and their literal root-deletion ports.

This inequality applies to the selected narrow facet. It does **not** constrain the old whole-polytope dual frequencies. Turning it into a dimension exclusion requires a facet-based weighting argument with a marked facet root, or another link between the distinguished edges and the small two-faces. Neither follows merely from the existence of a narrow facet.

Generation and saved replay passed. The source, certificate and report are frozen.

```sh
python3 scripts/fresh_sqrt5_narrow_e8_ports_20261016.py
python3 scripts/fresh_sqrt5_narrow_e8_ports_20261016.py --verify
```
