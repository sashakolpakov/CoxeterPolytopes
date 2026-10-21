# The remaining-seven contraction audit and the Möbius-plane inverse core

**Status (2026-08-09): exact seven-row audit and one new directional
upgrade.**  None of the seven inputs is a nonzero extension of the two
four-point `CI4` cores.  One input instead has a different, intrinsically
recognized finite-geometry core: the Möbius plane

\[
                         S(3,4,10)=\mathbf P^1(\mathbf F_9)
\]

with its thirty `F3`-sublines.  The general Frobenius-graph construction
realizes this core, the divisor formulas single out `q=3`, and the existing
uniform Coble--Vinberg theorem returns the exact same 206-wall polytope.

This upgrades that calibration from `SURFACE -> POLYTOPE` to a structural

```text
POLYTOPE -> finite geometry -> strong Coble surface -> the same P
```

result.  It is **L3**, not formal `L4`: this extractor itself audits the
seven archived raw inputs.  A separate normalized arbitrary-input recognizer
now covers the five-core union, but its rigid grammar gives no nonthin
pushback.

The exact artifacts are:

- [certificate](../scripts/certify_coble_remaining_seven_contraction_designs.py);
- [JSON](../data/enumerations/coble_remaining_seven_contraction_designs.json);
- [tests](../tests/test_coble_remaining_seven_contraction_designs.py).

**Post-audit promotion.**  The 78-wall Desargues candidate isolated below is
now closed independently.  The characteristic-five quotient plane realizes
all ten triangle lines and ten cut conics; cubic evaluation is injective; the
unique sextic double at all ten centres is geometrically integral with
exactly ten ordinary nodes; and the determinant-one raw marking returns the
same 78-wall chamber.  See the
[Desargues inverse theorem](coble_type_vi_char5_desargues_inverse.md).  The
seven-row JSON remains the frozen discovery audit, so its historical
`not_promoted` status describes this artifact at the time it was generated.

## 1. Input policy

For each row the extractor reads only

```text
integral lattice Gram;
primitive wall vectors.
```

It identifies boundaries by norm and characteristic parity, sets
`K=-sum(B_i)/2`, partitions the other walls by square and boundary profile,
and performs contractions using

\[
 D'.F'=D.F+(D.E)(F.E)\qquad(E^2=-1).                 \tag{1.1}
\]

No surface name, curve label, stored contraction, or known incidence
configuration participates in recognition.  Archive names occur only in the
human-readable lineage comparison after extraction.

## 2. Why none of the seven is a `CI4` extension

The obstruction is invariant and immediate once stated correctly.

- The complete conic polar/chord core has fourteen blowup centres and hence
  Picard rank at least `15`.
- The complete graph-projectivity core has sixteen blowup centres on
  `P1 x P1` and hence Picard rank at least `18`.
- An extension adds centres or curves; it cannot lower Picard rank.

Every remaining input has Picard rank `11` or `12`.  Therefore none can
contain either complete core, with or without a nonzero extension.  This is
not a failure to find a preferred contraction: it excludes every contraction
branch by rank.

Consequently `CI4` remains only the exact two-member core calibration and is
rated `L3`.  In particular, the former broad “standard-free extension” lemma
must not be used: base-point-free and generically separable do not by
themselves exclude a system composed with a pencil.  A future extension
verifier must require either

1. an explicit geometrically integral section certificate, or
2. after resolving the assigned base cluster, image dimension at least two,
   separability, and the remaining correct Bertini/avoidance hypotheses.

No nonzero extension satisfying that repaired gate is presently certified.

## 3. Exhaustive seven-row outcome

The simultaneous column records the largest set of pairwise disjoint
`(-1)` walls visible before any contraction.  For the first four rows the
script's optional full mode enumerates every reachable rank-one quotient,
not merely greedy branches.

| raw numerical row | rank / walls | largest simultaneous set | exhaustive `P2` targets | endpoint degree signatures | outcome |
|---|---:|---:|---:|---:|---|
| one boundary, 22 walls | `11/22` | 4 | 74 | 16 | sequential and noncanonical |
| one boundary, 30 walls | `11/30` | 7 | 573 | 28 | sequential and noncanonical |
| two boundaries, 26 walls | `12/26` | 4 | 52 | 9 | sequential and noncanonical |
| two boundaries, 36 walls | `12/36` | 7 | 890 | 28 | sequential and noncanonical |
| one boundary, 78 walls | `11/78` | 10 | one simultaneous plane frame | one | Desargues/cut design; explicit carrier completed in the post-audit theorem |
| one boundary, 72 walls | `11/72` | 10 | five simultaneous plane frames | one | tetrahedral `K6`-factorization candidate |
| two boundaries, 206 walls | `12/206` | 11 | 154 plane frames, but one intrinsic ruled frame | one canonical ruled core | Möbius-plane theorem |

The multiplicity of sequential endpoint signatures in the first four rows
freezes a useful obstruction taxonomy:

```text
RANK:       too small for either CI4 core;
FRAME:      no full simultaneous exceptional frame;
BRANCH:     many sequential endpoints survive;
DESIGN:     their boundary/root degree data have 9--28 distinct signatures.
```

This is not a surface nonexistence theorem.  It says that those rows do not
provide a canonical finite projective design under the present numerical
extractor, so no uniform design-realization theorem is claimed for them.

## 4. Label-free extraction of the 206-wall core

Start with its raw rank-twelve Gram and 206 primitive walls.  Norm and
boundary-pairing tests give

```text
2  square-minus-four boundaries;
30 square-minus-two walls;
10 square-minus-one walls with boundary profile (1,1);
164 square-minus-one walls with boundary profile (0,2).
```

The ten `(1,1)`-profile walls are pairwise orthogonal.  Contract them.  If
`D_0,D_1` are the pullbacks of the two boundary images, their Gram is

\[
                       \begin{pmatrix}6&10\\10&6\end{pmatrix}.
\]

The two classes

\[
             F={3D_1-D_0\over8},\qquad
             G={3D_0-D_1\over8}                         \tag{4.1}
\]

are integral and satisfy

\[
 F^2=G^2=0,\qquad F.G=1.
\]

Together with the ten exceptional walls they form a determinant-one basis
of the raw lattice.  Thus the endpoint is intrinsically `P1 x P1`; no
surface model was supplied.

Every square-minus-two wall is exactly

\[
                   F+G-\sum_{i\in A}E_i                \tag{4.2}
\]

for one four-subset `A` of the ten centres.  The thirty recovered blocks
satisfy

```text
each point occurs in 12 blocks;
each pair occurs in 4 blocks;
each triple occurs in exactly 1 block.
```

Hence the wall lattice itself produces a `3-(10,4,1)` Steiner system.  An
exact colored-incidence isomorphism identifies it with the Möbius/inversive
plane of order three on `P1(F9)`.  The particular isomorphism is deliberately
not archived, because the 1,440 design automorphisms make such a map
noncanonical.

## 5. The general Möbius-plane factory and its rigidity

Let `q` be a prime power and work over an algebraically closed field
containing `F_(q^2)`.  On

\[
                         Q=\mathbf P^1\times\mathbf P^1
\]

take the graph of `q`-Frobenius and its transpose:

\[
 B_0:\ y=x^q,\qquad B_1:\ x=y^q.
\]

They are smooth rational curves of bidegrees `(1,q)` and `(q,1)`.  Their
intersection is transverse and consists of

\[
                         q^2+1
\]

points indexed by `P1(F_(q^2))`.

The images of `P1(F_q)` under `PGL_2(F_(q^2))` form the Möbius design

\[
                 3-(q^2+1,q+1,1),
\]

with `q(q^2+1)` circles.  The twisted action

\[
               (x,y)\longmapsto (g(x),g^{(q)}(y))
\]

preserves the Frobenius graph and carries the diagonal `(1,1)` curve through
`P1(F_q)` to a unique smooth `(1,1)` curve through every circle block.

Blow up the `q^2+1` intersection points and write `E` for their sum.  The
strict transforms are

\[
 B_0=F+qG-E,\qquad B_1=qF+G-E,
\]

while a circle has class

\[
 C_A=F+G-\sum_{a\in A}E_a.
\]

Direct intersection theory gives the four decisive formulas

\[
\begin{aligned}
 B_i^2&=-(q-1)^2,\\
 C_A^2&=1-q,\\
 K.C_A&=q-3,\\
 B_0+B_1&=-2K+(q-3)(F+G).                \tag{5.1}
\end{aligned}
\]

Therefore

\[
\boxed{
 \begin{array}{c}
 B_0,B_1\text{ are disjoint }(-4)\text{-boundaries summing to }-2K,\\
 C_A\text{ are adjunction-compatible }(-2)\text{-roots}
 \end{array}
 \iff q=3.}
                                                               \tag{5.2}
\]

This is the theoretical gain over a named-row reconstruction: an infinite
finite-geometry factory has an exact necessary-and-sufficient Coble/root
rigidity criterion, and the raw polytope determines its accepted `q=3`
member.

## 6. Strong-Coble construction at `q=3`

At `q=3`, blow up the ten points of `P1(F9)` on the two Frobenius graphs.
Equation (5.1) gives two disjoint `(-4)` curves with sum `-2K` and thirty
actual `(-2)` curves (4.2).

The anti-canonical system is empty.  A section of `-K` would be a `(2,2)`
curve through all ten points.  Restriction to the Frobenius graph identifies
the nine bidegree monomials with all monomials of degrees `0,...,8` on
`P1`; the exact ten-point evaluation has rank nine and zero kernel.

The anti-bicanonical boundary is unique.  A `(4,4)` divisor double at all ten
centres has intersection 16 with either Frobenius graph, but local
intersection at the ten double points is at least 20.  It must contain that
graph.  The residual `(3,1)` divisor passes through all ten points; since its
intersection with the other `(3,1)` graph is only six, it must equal that
graph.  Thus

\[
                  |-K|=\varnothing,\qquad
                  |-2K|=\{B_0+B_1\}.
\]

The resulting rational surface is a strong Coble surface constructed from
the P-derived incidence core, without importing classified surface
existence.

## 7. Exact return to the same 206-wall polytope

The determinant-one marking identifies the constructed boundaries, ten
crossing curves, and thirty roots literally with the raw classes.  The input
already passes the surface-free `D_CobVin` gates.  The general theorem then
applies:

1. Vinberg finite index and Kondō's root-exhaustion theorem exclude hidden
   Coble--Mukai roots.
2. The primitive-isotropic/Halphen converse produces exactly the 164
   integral `(0,2)` cusp lifts and excludes hidden `(-1)` curves.
3. Therefore all negative curves and all nef inequalities are exactly the
   206 raw walls.

Hence

\[
                 \boxed{\mathbf P(\operatorname{Nef}S)=P}
\]

under the original integral marking and orientation.  Finite volume remains
a P-side fact from the already certified exact cone/Vinberg audit; it is not
inferred from finite sidedness.

Define `D_CobVin^Mob3` structurally as the members of `D_CobVin` whose
intrinsic two-boundary crossing contraction is `U` and whose root-centre
hypergraph is the Möbius `3-(10,4,1)` design.  The preceding argument proves
nonemptiness and exact return for every input accepted by an
arbitrary-input implementation of that recognizer.  The current executable
certifies the archived calibration but does not yet expose that normalized
interface, so the production rating remains `L3`.

## 8. The other two finite designs

The audit also recovers two additional designs.  Both branches have since
been completed by explicit P-first constructions.  Neither is promoted merely
by importing a known surface.

### The 78-wall plane design

There is one simultaneous plane frame.  Its ten centres identify with the
edges of `K5`.  The roots push down to

```text
10 lines through the three edges of each K5 triangle;
10 conics through the six-edge cut delta(A), A a two-subset of K5.
```

The line part is the Desargues `10_3` configuration.  The formerly missing
certificate is now frozen: an explicit geometrically integral ten-node
sextic, determinant-one raw basis, twenty vector-level root identities, and
120 colored incidence isomorphisms give an exact P-first same-chamber return.
See [`coble_type_vi_char5_desargues_inverse.md`](coble_type_vi_char5_desargues_inverse.md).

### The 72-wall plane design

There are five plane frames, all with the same profile.  Each splits the ten
centres into four special points and six conic points.  The twelve lines give
four perfect matchings of `K6`; one conic contains all six ordinary points,
and three more use the pairs of the fifth matching together with all four
special points.  Four quartics are each double at three special points.

This candidate is now exact over `F5`.  The six ordinary points lie on
`xz-y^2=0`; four concurrent perfect matchings give the four special poles.
The twelve lines, four smooth conics, and four geometrically integral nodal
quartics realize all twenty roots.  The unique boundary sextic is
geometrically integral with four `A1` nodes and six `A2` cusps--not ten
nodes--and its strict transform is smooth.  All five intrinsic frames and all
120 colored markings return the exact root vectors, the twenty-root cone has
`450+51+0` rays with its 51 cusps exactly `{E+b}`, and all 72 walls return.
See [`coble_type_vii_char5_tetrahedral_inverse.md`](coble_type_vii_char5_tetrahedral_inverse.md).

## 9. Lineage

Kondō's *Coble surfaces with finite automorphism group*, Example 4.17,
already publishes the characteristic-three construction on `P1 x P1`, the
two `(1,3)/(3,1)` boundaries, their ten `F9` intersection points, and the
thirty `(1,1)` curves.  Kondō credits the construction to Katsura--Kondō,
Section 4, and says the surface was suggested by Mukai and Ohashi.

Accordingly, this report does **not** claim that the `F9` surface, the thirty
curves, or the abstract 206-wall polytope is new.  The project-new scope is:

- raw-P, label-free recovery of the ruled contraction;
- intrinsic recognition of the thirty blocks as `S(3,4,10)`;
- the general `q`-factory rigidity statement (5.2);
- construction followed by the exact same-P compiler, upgrading
  directionality.

## 10. Exact next gate and replay

The finite promotion gate was a normalized input schema

```text
{ lattice_gram, primitive_wall_vectors, orientation/interior_vector }
```

and one arbitrary-input recognizer for all five finite-geometry cores:

```text
polar/chord CI4;
grid CI4;
Mobius 3-(q^2+1,q+1,1).
K5 Desargues triangle/cut design.
```

That five-core interface is now frozen in
[`coble_finite_geometry_raw_recognizer.md`](coble_finite_geometry_raw_recognizer.md):
it neither selects archive paths nor compares an output Gram with a stored
answer.  It closes generic recognition but not the theorem payoff.  The
rigid five-core list gives only thin Kondō pushback.  The real promotion gate
is now a
nonzero P-only extension class with uniformly nonempty incidence/Hilbert
moduli and a classification consequence not already encoded by the core.

Quick replay:

```bash
python3 scripts/certify_coble_remaining_seven_contraction_designs.py --check-archive
python3 -m pytest -q tests/test_coble_remaining_seven_contraction_designs.py
```

The four exhaustive sequential censuses are optional because the 36-wall
replay is substantially slower:

```bash
python3 scripts/certify_coble_remaining_seven_contraction_designs.py \
  --full-contractions --check-archive
```
