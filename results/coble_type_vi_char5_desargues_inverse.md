# A raw-polytope Desargues construction of the characteristic-five `H10/78` chamber

**Status (2026-08-09): exact structural `POLYTOPE -> SURFACE -> same P`
calibration.  Inverse level `L3`, not `L4`.**

Starting only from the rank-eleven lattice and 78 primitive walls of the
finite-volume chamber `COBLE-VI-p5-b1-H10-78`, the contraction extractor finds
a unique simultaneous set of ten disjoint `(-1)` walls.  Contracting them
gives `P2`.  The twenty square-minus-two walls become

```text
10 lines through triples of centres;
10 conics through six centres.
```

Their incidence is intrinsic: the centres are the ten edges of `K5`, the
line triples are its ten triangles, and the conic blocks are its ten six-edge
cuts.  This note realizes that P-derived design in characteristic five,
constructs its unique anti-bicanonical boundary, proves that the result is a
strong Coble surface, and invokes the already frozen root/cusp theorem to
return literally the original 78-wall chamber.

The result upgrades this row from a surface-first realization to a genuine
P-first carrier theorem.  It does **not** produce a new abstract Coxeter
polytope, exclude alternative carriers, or prove a class-wide dimension
bound.

The exact artifacts are:

- [certifier](../scripts/certify_coble_type_vi_char5_desargues_inverse.py);
- [JSON](../data/enumerations/coble_type_vi_char5_desargues_inverse.json);
- [tests](../tests/test_coble_type_vi_char5_desargues_inverse.py);
- [raw seven-row extractor](coble_remaining_seven_contraction_design_audit.md).

## 1. The characteristic-five quotient plane

Over `F5`, put

\[
 W=\{(x_0,\ldots,x_4):\sum x_i=0\},\qquad
 V=W/\langle(1,1,1,1,1)\rangle.
\]

The all-ones vector lies in `W` precisely because the characteristic is five,
so `dim(V)=3`.  For every unordered pair `{i,j}`, let

\[
                         p_{ij}=[e_i-e_j]\in\mathbf P(V).
\]

The ten projective points are distinct.  For every triple `{i,j,k}`, the
identity

\[
              (e_i-e_j)+(e_j-e_k)+(e_k-e_i)=0
\]

puts `p_ij,p_jk,p_ik` on one line.  Exact finite-field row reduction gives
one and only one line through each of the ten triangle triples.  Likewise,
for each two-subset `A` of the five vertices, the six points indexed by the
cut `delta(A)` impose rank five on the six-dimensional conic system.  Hence
they determine one conic.  All ten conics have nonsingular symmetric matrices.

The resulting root classes on the blowup are

\[
 L_T=H-\sum_{e\in T}E_e,
 \qquad
 Q_A=2H-\sum_{e\in\delta(A)}E_e.
\]

Thus

\[
 L_T^2=Q_A^2=-2,
 \qquad K.L_T=K.Q_A=0.
\]

Their complete degree/multiplicity profile is exactly the one independently
recovered from the raw wall lattice.

## 2. The unique sextic boundary

Let `Z` be the ten Desargues points.  Cubic evaluation has full rank

```text
dim H0(O_P2(3)) = 10,
rank(ev_Z)      = 10,
```

so no cubic passes through every point and `|-K|` is empty after blowing up
`Z`.

For sextics, the thirty first-derivative conditions for a double point at
every member of `Z` have exact rank `27` in the 28-dimensional space of
degree-six forms.  Their kernel is one-dimensional.  In the certified
coordinates its normalized generator is

\[
\begin{aligned}
F={}&x^4y^2-x^4yz+x^4z^2+2x^3y^3+2x^3z^3\\
 &+x^2y^4+x^2y^2z^2+x^2z^4-xy^4z-xyz^4\\
 &+y^4z^2+2y^3z^3+y^2z^4
 \quad\in\mathbf F_5[x,y,z].
\end{aligned}
\]

Exact Singular factorization has only the unit and `F`, with exponent one;
therefore `F` is irreducible over `F5`.  Saturating its Jacobian ideal by the
irrelevant ideal gives a projective zero-dimensional scheme of degree ten.
The ten prescribed `F5` points are already distinct points of that scheme,
so they are its complete reduced support.  Direct quadratic-tangent-cone
calculation has nonzero discriminant at all ten points.  Hence `F=0` has
exactly ten ordinary nodes and no other geometric singularities.

Irreducibility over `F5` also implies geometric irreducibility here.  If `F`
split over the algebraic closure into `r>=2` Frobenius-conjugate components,
each `F5`-rational node would lie on every component.  Their common degree
would be `6/r<=3`, so two distinct components would share ten points,
contradicting Bézout because `(6/r)^2<=9`.

Finally, all 90 incidences of a triangle line or cut conic with a prescribed
centre have tangent direction different from both branches of the node.
Bézout then exhausts the intersections:

```text
line . sextic  = 6  = 3 centres * multiplicity 2;
conic . sextic = 12 = 6 centres * multiplicity 2.
```

After one blowup at each centre, the boundary is therefore smooth and
disjoint from all twenty roots.

## 3. Strong Coble surface and exact return

On `S=Bl_Z(P2)`, the strict sextic has class

\[
              B=6H-2\sum E_i=-2K_S,
              \qquad B^2=36-40=-4.
\]

It is smooth, geometrically integral, and rational: a plane sextic has
arithmetic genus ten, exactly cancelled by the ten nodes.  The previous
linear-system calculations give

\[
                         |-K_S|=\varnothing,
              \qquad |-2K_S|=\{B\}.
\]

Thus `S` is a strong Coble surface.  More importantly, the construction is
marked by the raw polytope: its unique ten-wall contraction frame is the
exceptional basis, and its twenty root supports are precisely the triangle
and cut blocks realized above.

The marking is checked at vector level, not inferred from the degree
histogram.  In raw wall order the unique exceptional frame is

```text
25, 26, 29, 35, 40, 43, 49, 57, 64, 77.
```

If `b=-K` is the unique characteristic norm-minus-one wall, then
`H=(b+sum E_i)/3` is integral and

```text
det(H,E_1,...,E_10) = +/-1,
Gram(H,E_1,...,E_10)=diag(1,-1^10).
```

Every one of the twenty raw root vectors is exactly
`dH-sum_(i in A)E_i`, not just Gram-isomorphic to such a class.  The complete
colored point/line/conic incidence graph has exactly 120 isomorphisms to the
`K5` triangle/cut design; a lexicographically canonical one is frozen in the
archive.  This supplies the exact marked root-seed transport required below.

The general `D_CobVin` return theorem now applies without importing surface
existence.  Finite-index reflection theory exhausts the twenty effective
Coble--Mukai roots; the 57 primitive cusps lift to all 57 `(-1)` curves; and
the boundary, roots, and cusp lifts recover every original wall.  Hence

\[
                 \mathbf P(\operatorname{Nef}S)
                    =P_{\mathrm{VI},p=5}
\]

as a marked chamber, with the already certified exact ray census

```text
21,661 finite + 1,797 ideal + 0 hyperideal.
```

The literal Vinberg edge and parabolic-neighbourhood audits remain those of
the frozen full-nef theorem; no finite-volume inference is made merely from
the construction of finitely many curves.

## 4. Scope and lineage

Kondō's Type-VI characteristic-five Coble surface, its twenty
`Gamma_VI` roots, `S5` symmetry, and uniqueness are known.  The complete
78-wall full-nef chamber was already a project-positive result.  The new
content here is directional and constructive:

```text
raw 78-wall P
  -> unique P2 contraction
  -> intrinsic K5 triangle/cut design
  -> explicit integral sextic Coble carrier
  -> exact same 78-wall P.
```

Worldwide priority for this plane presentation or the literal subgroup
chamber remains open, so neither is claimed.  The inverse rating is `L3`:
one structural P-first calibration now closes, but the result does not prove
that an independently broad class of polytopes has exhaustive surface
carriers, nor does independent AG yet push a new bound back to that class.

Primary lineage: Shigeyuki Kondō,
[*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591),
especially Theorem 1.1, Propositions 3.1, 3.7, 3.8, Lemma 3.11, and
Example 4.12.

## 5. Replay

```bash
python3 scripts/certify_coble_type_vi_char5_desargues_inverse.py --check-archive
python3 -m pytest -q tests/test_coble_type_vi_char5_desargues_inverse.py
```
