The largest admissible exceptional-cap heights give a changed H9 seed with **15 actual facets, exactly62 finite vertices and220 finite edges**. The richer finite contacts do not close the chamber: an exact strict spacelike escape proves infinite volume. Moreover, the same escape proves infinite volume for **every independent positive height choice on these three fixed projected cap directions**, retaining the source walls and bottom.

The [producer](../scripts/fresh_h9_h8_cap_height_20261014.py) and [exact certificate](fresh_h9_h8_cap_height_20261014.json) use the frozen [compact H8 source](../data/known_examples/bugaenko_compact_h8.json) and the independently reviewed [earlier15-wall lift](fresh_h9_h8_disjoint_cap_20261014.md). No frozen source or earlier construction is changed.

Write \(a=(\sqrt5-1)/4\), \(\phi=(1+\sqrt5)/2\). The three selected future source rays \(v_0,v_{32},v_{33}\) have square \(-a\) and minimum nonzero source-wall contact1. Their nonzero contact lists are

| Source vertex | Nonzero contacts \(-(v_i,u_j)\) |
|---|---|
|0|\(\phi,1,1\)|
|32|\(\phi,1,\phi^4\)|
|33|\(1,\phi,\phi^4\)|

For a general positive vertical coefficient \(d\) and cap \((v_i,-h_i)\), \(h_i>0\), let \(q_i=dh_i^2-a>0\) be its squared norm. The source contact1 has normalized cosine \(1/\sqrt{2q_i}\). Therefore its complete Coxeter-compatible possibilities are

\[
0<q_i\le\tfrac12,
\qquad\text{or}\qquad
q_i=\frac1{2\cos^2(\pi/m)},\quad m\ge3.
\]

In particular \(q_i\le2\), with no label cutoff. Conversely these possibilities pass all source-wall gates. For \(m=3\), the contact \(\phi\) becomes \(\pi/5\); for every \(m\ge4\), it is strictly dotted since \(\phi/\sqrt2>1\). The contact \(\phi^4\) is always dotted. At \(q_i=1/2\), minimum contacts are parallel; below it they are dotted.

The bottom remains strictly dotted to every cap, because its pairing discriminant is \(d a>0\). The three caps are mutually strictly dotted for every allowed independent height. Indeed, writing \(b_{ij}=-(v_i,v_j)\),

\[
(S_i,S_j)=-b_{ij}+\sqrt{(a+q_i)(a+q_j)}.
\]

Each exact \(b_{ij}>a+4\), while
\(\sqrt{(a+q_i)(a+q_j)}+\sqrt{q_iq_j}\le a+4\).
Thus their mutual pairing is less than \(-\sqrt{q_iq_j}\). Independent heights cannot create finite contacts between the three caps.

For the maximal choice \(q_i=2\), use

\[
D=a+2=(7+\sqrt5)/4,\qquad Q=Q_8\oplus\langle D\rangle,
\qquad S_i=(v_i,-1).
\]

The full105-pair census is79 right,12 at \(\pi/3\),5 at \(\pi/5\), and9 dotted. An exact strict center and its orthogonal projections certify all15 actual facets. The41 original bottom vertices survive. The three caps respectively have11,5 and5 finite vertices, giving62 in total. All123 combinations of a cap and a source vertex are checked exactly, including feasibility and causal type; none has zero norm.

This finite-vertex list is complete. A finite vertex can contain at most one of the bottom and the three mutually dotted caps. A source-only set cannot give nine positive independent normals in the source Lorentz space. Projection of a finite vertex containing one cap to the source is timelike and must lie on eight source walls, hence is one of the41 source vertices. Each actual finite vertex has exactly nine incident facets. Their complete shared-eight-facet list gives220 finite edges and \(62\cdot9-2\cdot220=118\) remaining finite-corner edge germs. The latter are not claimed to have bounded ideal endpoints.

The uniform escape has a simple literal formula, using the archived source center \(H_8\):

\[
E=(v_0+v_{32}+v_{33}+H_8/16,-5),\qquad
E^2=\frac{4363}{256}-\frac{3057\sqrt5}{512}>0.
\]

All15 wall pairings are strictly negative, as is its pairing with the strict future center. Thus the maximal chamber has an open ideal-boundary region.

For arbitrary \(d>0\) and positive heights, an isometry rescales the vertical coordinate to this same form \(Q\). The caps become \((v_i,-\eta_i)\), where

\[
\eta_i=\sqrt{\frac{a+q_i}{D}}\le1.
\]

Their pairings with this same escape satisfy

\[
(E,(v_i,-\eta_i))-(E,(v_i,-1))=5D(\eta_i-1)\le0.
\]

Consequently every admissible independent height choice, and every subfamily of these caps, retains the strict escape. The argument is uniform over real heights and all finite Coxeter labels. It also shows why increasing the finite-vertex count at the maximal height does not produce finite volume.

The maximal coefficient has field norm11/4. At the other embedding of \(\mathbb Q(\sqrt5)\), the form is positive definite. These checks neither identify a reflective lattice nor prove nonreflectivity; no unimodular-host theorem is extrapolated to this form. The infinite-volume conclusion comes from the explicit escape.

Generation and saved replay pass. The theorem fixes the eleven source directions, the bottom, and these three projected cap rays. New projected cap directions, coupled changes to the source walls, and additional facets remain outside it.
