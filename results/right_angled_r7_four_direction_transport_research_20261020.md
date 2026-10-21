# A second actual transport: four central R7 faces

20 October2026. **New paper mechanism for the mixed dimension-ten problem.** No scientific job, census, optimizer, checker, Lean or build ran. This proves an actual reversible move and ordinary defect/finite-edge budgets, not a new dimension ceiling. The existing dimension-ten theorem and its manuscript remain unchanged.

Use the actual P3, Q8, R7 and small-three-face statements in `manuscripts/right_angled_dimension10.tex`. A finite vertex q in dimension10 has q²=-1 and orthonormal outward normals n1,...,n10. Direction i omits ni. In R7 the four finite vertices form a three-leaf star: the center has no ideal-ending edges, each leaf has two, and the three cusps form the ideal triangle. All three polygons at its center have three finite vertices and one cusp.

## 1. The exact central R7 metric

Take its center q and three directions i,j,k. Write Li for the finite edge lengths. The other endpoint and endpoint wall of direction i have the actual formulas

    qi=cosh(Li)q-sinh(Li)ni,
    ri=sinh(Li)q-cosh(Li)ni.

In the quadrilateral at directions i,j, the other two edge walls ri,rj are asymptotic at its ideal corner, so their unit product is-1. Since ni,nj are orthogonal, `sinh(Li)sinh(Lj)=1`. The three pair identities force every sinh(Li)=1. Thus

    qi=sqrt2 q-ni,       ri=q-sqrt2 ni,
    uij=sqrt2 q-ni-nj,
    tijk=sqrt2 q-ni-nj-nk.                       (1)

The uij are the three actual cusps. The last normal is the actual opposite ideal-triangle facet: it is unit, pairs to zero with all three cusps, and its product with q is negative. Its sign is fixed. Together with the three ni and three ri, it is the complete seven-wall R7 system. This derivation uses an actual face with its known incidence, not a hypothetical metric assigned to a count pattern.

## 2. Four central triples give an actual finite opposite

Let A be four finite-ending directions at q, and suppose all four of its three-subsets are actual R7 faces with q their center. Equation(1) gives all four actual walls

    ti=sqrt2 q-sum_(j in A\{i}) nj,    i in A,

and all ri=q-sqrt2 ni. The ti are mutually orthonormal and perpendicular to every old nl outside A. Define

    qA=3q-sqrt2 sum_(i in A) ni.                 (2)

It is future unit timelike, has product zero with all ti and outside nl, and product-sqrt2 with each removed ni.

Here is the full-wall feasibility argument. Any other actual wall not incident to q has the complete-frame form `r=a q-sum b_l nl`, with a>0,b_l>=0. If it is distinct from the four ti, acuteness with each ti gives

    sum_(j in A\{i}) bj <= sqrt2 a.

Summing gives `3 sum_A bj<=4sqrt2 a`, and therefore

    r.qA=-3a+sqrt2 sum_A bj<=-a/3<0.             (3)

The four ti and old incident normals have already been checked separately. Thus every wall is feasible, with precisely the ten displayed incident walls. In particular qA is an **actual simple finite vertex**, not a formal Lorentz point. Different four-sets give different qA by independence of the old normals.

## 3. Reversal and the actual R7 three-faces at qA

For each i in A the old actual vertex qi has bad triangular polygons on each pair of the other three A directions, by its three R7 leaf incidences. The actual P3 lemma therefore supplies

    q'i=2qi-sum_(j in A\{i}) nj
       =sqrt2 qA-ti.                            (4)

Its incident normals are ri, the three tj with j!=i, and old normals outside A. Hence qA and q'i share exactly nine incident walls and are endpoints of an actual finite edge. Their segment is feasible by convexity; no additional wall can create an intervening vertex while strict at both endpoints.

The reverse finite-edge endpoint walls and cusps are

    qA-sqrt2 ti=ri,
    sqrt2 qA-ti-tj=uij.                         (5)

For a reverse triple A\{l}, its opposite ideal-triangle wall is

    sqrt2 qA-sum_(i in A\{l}) ti=nl.            (6)

All seven reverse walls are therefore actual old walls, and its seven R7 vertices are actual points supplied by(2),(4),(5). In its Lorentz four-space these seven walls are exactly the standard R7 system(1). Their convex hull is finite-volume and all its vertices are feasible, so another wall cannot cut it. Thus the whole reverse three-face is R7, not merely a corner with the same ideal degree. Every reverse triple is central R7. Applying the construction again returns q.

We therefore have an undirected graph of actual finite vertices, with one edge for each such good four-set. This opposition is a diagonal of an actual four-face, not an ordinary polytope edge. No equality between this graph and the earlier P3 opposition graph is asserted.

## 4. Invariant, ratios and reciprocal transport

The sum of the four new normals in(2) proves that

    H=2sqrt2 q-sum_(i=1)^10 ni,     H²=2         (7)

is invariant along this graph. Equivalently the norm2 vector `vA=sqrt2 q-sum_A ni` defines the reflection `x -> x-(x.vA)vA`, sending the full frame at q to that at qA and fixing H. This is a frame isometry; it is not claimed to be an ambient polytope reflection or an arithmetic structure.

On each finite connected component choose one future null Z perpendicular to H, and put F=-Z.q>0,p_i=-Z.ni/F. Then

    sum p_i=2sqrt2,    sum p_i²=1,
    rA=F(qA)/F(q)=3-sqrt2 sum_A p_i.            (8)

The same Z is used at both ends; reverse ratios are reciprocal. Centering at sqrt2/5 gives the uniform formal range

    (7-2sqrt6)/5 <= rA <= (7+2sqrt6)/5,
    2/5 < rA < 5/2.                            (9)

For the actual moves, summing the endpoint terms `rA^-1-1` over all vertices is nonnegative: each undirected edge contributes `r+r^-1-2>=0`. This is an ordinary global transport statement. It neither averages pointwise maxima nor substitutes different null vectors at opposite endpoints.

For orientation, at the full ten-direction level the formal ratio has mean7/5 and variance8/75. The exact reciprocal identity

    1/(mu+d)=1/mu-d/mu²+d²/(mu²(mu+d))

and(9) give the unconditional rational bound

    average r^-1 <=5/7+(8/75)/((2/5)(7/5)²)
                 =125/147.

Hence the complete210 formal four-subsets satisfy

    sum_(|A|=4)(rA^-1-1) <=-220/7.              (10)

Deleting one unavailable move can increase this sum by at most3/5, since every term is greater than-3/5. This is already a negative row for a zero-ideal-degree vertex, with no new finite positivity certificate. Bounds for other finite-direction subset sizes are a separate analytic task; equation(10) alone does not close the global count.

## 5. Ordinary missing-move budget

Let k(q) be the number of ideal-ending directions and m(q)=10-k(q). Let z(q) count three-faces through q whose three directions are finite-ending but which are not central R7. The small-three-face classification shows every such incidence lies in a non-P3 three-face with nu>=11. For each four-set of finite directions that is not good, at least one of its four triples is counted by z(q). Thus, writing M(q) for the number of actual R7 four-moves,

    0<=binom(m(q),4)-M(q)<=(m(q)-3)_+ z(q).      (11)

This is a union bound on actual direction sets, with no classification of larger four-faces.

Let Z0=sum z(q). For a non-small three-face Y write e=nu(Y)-10>=1 and d=A2-3A0. The already proved one-cusp lemma gives p<=10e at nu11; for nu>=12 it follows from p<=2nu-4. Also

    4A0=p-d-A1-A3<=10e-d.

The three small types have d0, so summation gives

    4Z0<=10E-D,
    sum_q[binom(m(q),4)-M(q)]<=7Z0.            (12)

For any independently proved subset upper Um for the full formal sum in(8), equations(9),(11) and ordinary reciprocity yield the necessary inequality

    0<=sum_q Um(q)+(3/5)sum_q(m(q)-3)_+z(q)
      <=sum_q Um(q)+(21/20)(10E-D).            (13)

One may instead retain the separate z(q) distribution for a sharper count. A common rough error bound must not be mistaken for proof that the entire analytic row is strong enough to exclude dimension10.

## 6. A second ordinary constraint: regular finite edges flip degrees

Call an actual finite--finite edge regular if **every** three-face containing it is one of P3,Q8,R7. At one endpoint q exclude its edge direction and consider any pair of the other nine directions. If both are finite-ending, its three-face is central R7, and both directions become ideal-ending at the other endpoint. If both are ideal-ending, it is a leaf of R7, and both become finite-ending. If one is of each kind, the face is Q8 and exactly one remains ideal-ending at the other endpoint.

Among nine transverse directions, one class has at least two members. The first or second observation fixes that class; the mixed observation then fixes all of the other class, even when it has size zero or one. Every transverse status is therefore reversed. Consequently

    k(q')=9-k(q).                              (14)

Such an edge also has sinh(length)=1: some same-status transverse pair gives an actual R7 containing it. This metric statement is not required for the incidence count.

Let f_k count ambient finite vertices of ideal degree k, and let N_k count incidences at such vertices with irregular finite edges. Regular edges give exactly

    (10-k)f_k-N_k=(k+1)f_(9-k)-N_(9-k),
    0<=k<=9.                                  (15)

An irregular edge lies in a nu>=11 three-face. A union bound over those actual faces gives

    sum N_k <=sum_(big Y) 2e_FF(Y)<=30E.        (16)

More precisely `2e_FF=3p-e_FI=2p-d-2A0-2A3`, so

    sum N_k <=20E-D-2Z0.                       (17)

Both bounds are valid; one may use their minimum. In particular the sum over the five unordered pairs k,9-k of the absolute discrepancies in(15) is at most this same total, without an extra factor two. These are finite-edge constraints, independent of the two diagonal transport graphs and without any global no-finite--finite-edge assumption.

At E=0 equation(15) forces paired degree classes, including the formerly problematic exact ratio f9=10f0. It does not say those counts are geometrically realizable. The four-move transport now acts on the finite-direction structure of precisely this obstruction, while the old P3 transport acts on its ideal-direction structure. Their joint quantitative consequence is still to be proved; no LP, new scalar certificate, dimension-nine theorem or example has been obtained here.
