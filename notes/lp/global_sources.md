# Primary-source constraints for the global formulation

These reading notes retain the useful uniform statements from the discarded
search. They are inputs to the [global necessary-data formulation](global_formulation.md),
with a separate author review still required before accepting a new bound.
Source-byte bindings are in [sources.json](../../research/global_formulation/sources.json).
Third-party PDFs and extracted texts remain outside the repository.

## Sources and scope

- [Andreev, Math. USSR-Sb. 12 (1970), 255–259](https://www.mathnet.ru/eng/sm3511),
  supplied as `sm3511_eng.pdf`: all five pages read; Theorems 1–2 and the
  necessity argument checked against page images where OCR was ambiguous.
- [Prokhorov, Math. USSR-Izv. 28 (1987), 401–411](https://www.mathnet.ru/eng/im1484),
  supplied as `im1484_eng.pdf`: all eleven pages read. The local cases in
  Lemma 1, the displayed incidence formulas in Lemma 2, and Figures 3–4
  with Propositions 4–5 were inspected in page images.
- [Vinberg, Hyperbolic reflection groups (1985)](https://www.mathnet.ru/eng/rm2140):
  Theorems 3.1–3.2, Proposition 4.2 and Theorem 4.1 supply the face,
  cusp, endpoint and finite-volume compatibility background. The preceding
  audit used the English online text; a local full-PDF acquisition is not claimed.
- [Vinberg, The absence of crystallographic groups of reflections in
  Lobachevsky spaces of large dimension (1984)](https://www.mathnet.ru/eng/mmo444):
  the preceding audit read the Russian PDF. Proposition 8(5), p. 76,
  uses the quadrilateral-star determinant condition; Section 6 describes
  diagram colorings. Section 6.12, p. 91, gives the E7 warning about
  restricting the weight family, cited in the LP manuscript's remark.

The supplied PDFs were read in full in the preceding source audit.
The retained metadata records their byte hashes and the inspected scope.
This does not reprove the cited classifications or Khovanskii’s theorem.

## Andreev's finite-volume conditions

Theorem 2 excludes the simplex from its sufficiency statement. Its
necessity conditions, pp. 256–257, are:

| Condition | Content | Use in this program |
| --- | --- | --- |
| m0, m1* | Nonobtuse positive angles; vertex angle sums for valence three and four | Already largely in spherical/Euclidean vertex figures; retain the ideal equality cases. |
| m2, m3 | Strict angle-sum bounds on triangular and quadrangular prismatic elements | Apply with the actual incidence hypotheses, including ideal intersections. |
| m4 | In a triangular prism, the six base-to-side dihedral angles have total less than 3 pi | Excludes all six being right; useful as a shared three-face filter. |
| m5 | If first and third faces touch at a vertex, the middle face is adjacent to both, and there is no common triple vertex, the two intervening angles have total less than pi | A useful cusp-incidence filter beyond separate corner PSD tests. |

For an actual three-face of the original chamber, use its *projected and
normalized* Gram entries. Such a face is acute-angled and finite-volume,
but its intrinsic angles need not be Coxeter angles. Thus no discrete
angle margin may be assumed for projected three-face angles. For m4
and m5, nonobtuseness makes their strict sums equivalent to forbidding
the specified all-right pattern. An SDP using only non-strict weak
inequalities would lose that distinction. Retain incidence strata or
prove a quantitative margin before replacing strictness numerically.

Here is a direct version of m5 in every ambient dimension. Let actual
unit facet normals a,c satisfy <a,c>=-1 and let their facets touch at
the ideal point v. Their nonzero null sum a+c spans v. If another
facet normal b satisfies <a,b>=<b,c>=0, then <b,a+c>=0. Its hyperplane
also contains v; because v belongs to the closed chamber, v belongs
to all three closed facets. Consequently the incidence stratum with
no common triple cusp cannot have both intervening angles right.
The sum a+c is nonzero: otherwise the two opposite halfspaces would
leave the chamber with empty interior.

For *original Coxeter facets*, put u=-<a,b> and w=-<b,c>. Then u,w>=0
and each nonzero one is at least 1/2. The excluded right pattern has
the quantitative linear form

    u+w >= 1/2.

This uses the original Coxeter-angle discreteness and does not apply
with that margin after projection to an arbitrary three-face. The
matrix with diagonal 1 and off-diagonal entries g_ac=-1, g_ab=g_bc=0
is PSD of rank two. Matrix positivity by itself therefore misses the
required incidence consequence: the three facets must share the cusp.

## Prokhorov's local information before numerical aggregation

On p. 404, a bad three-face is exactly a triangular bipyramid. Its
three equatorial vertices are ideal; rectangular horospherical links
force every one of its nine intrinsic dihedral angles to be right.
Its two apices are finite and have right trihedra. This completes the
geometric identification left open in the earlier link note.

Write T=f_3 and B for all and bad three-faces. For a four-face with
k three-faces, let b be its number of bad ones. Lemma 1's actual
case proofs, p. 405, give the following caps, some stronger than its
single displayed estimate b<=2k-10:

| k | 5 | 6 | 7 | 8 | 9 | at least 10 |
| --- | --- | --- | --- | --- | --- | --- |
| cap on b | 0 | 0 | 3 | 6 | 7 | k |

Keep separate variables h_k and b_k for four-face counts and bad
incidences at k=5,...,9. All k>=10 can be retained through three
aggregate variables h_+, K_+, b_+, with K_+>=10h_+ and 0<=b_+<=K_+.
This is a finite outer representation of these selected conditions;
it imposes no maximum on k. Nonnegative aggregate directions with
h_+=0 are permitted by the closed relaxation, not claimed realizable.

With K_4=sum k h_k+K_+ and B_4=sum b_k+b_+, the exact marginals are

    K_4=(n-3)T,       B_4=(n-3)B,       f_4=sum h_k+h_+.

Each three-face has exactly n-3 containing four-faces: its relative
interior is a finite point and the local interval of containing faces
is Boolean. Nonsimple ideal vertices cause no extra factor here.

Prokhorov's Figure 3, pp. 409–410, is also useful for the Gram oracle.
The four shared marked-star types yield charges 3, 2, 1 and 1/3 for
his particular distance weight. Retain their common white block,
marked pairs, deletion types and overlapping paths. These are not
four new universally valid inequalities for an arbitrary replacement
weight vector. New weights require the corresponding shared-state
inequalities to be imposed on their actual feature sums.

## Coupling the exceptions to both vertex kinds

For a finite vertex v with spherical component ranks s_j define

    A_v=sum_{i<j<k} s_i s_j s_k.

Exactly A_v incident three-faces have a right trihedron at v. To see
this, release three nodes I of the spherical diagram and project away
the other nodes. The inverse of the residual Gram block is the
principal I-block of the full inverse Gram matrix. On a connected
spherical component that inverse has strictly positive entries, and
between components its entries vanish. The residual three normals
are mutually orthogonal exactly when the released nodes lie in three
different connected components. This proves the formula with actual
embedding multiplicity, irrespective of the Coxeter labels.

Both finite apices of a bad face consume different such incidences:

    2B <= A,           A=sum_{v finite} A_v.

For the diagnostic forest 2D_17+D_23, A_v=17*17*23=6647, whereas the
total number of incident three-faces is binomial(57,3)=29260. This is
a new statistic attached to the same full vertex column, not a new
independent choice for that column. This forest is still only formal
necessary data, not a realized example.

At an ideal vertex with simplex-factor ranks r_i, put R=n-1 and
V=product(r_i+1). The earlier quadrilateral-link count is

    Q_v=V*(R^2-sum r_i^2)/8,       3B <= Q=sum Q_v.

Let P count triangular link two-faces at *all* vertices. Its finite
column is binomial(n,3); its ideal column is V*sum r_i(r_i-1)/6.
Every vertex of a three-face has valence three or four, and its
vertex figure is the corresponding triangular or quadrilateral link
two-face. Double-counting edges and applying Euler to each three-face
give the exact further couplings

    2E_3=3P+4Q,         2K_3=4T+P+2Q,
    E_3=binomial(n-1,2)f_1,       K_3=(n-2)f_2.

These relations make the link capacities part of the incidence system.
P,Q,A must be computed from the *same* finite and ideal vertex counts
as the original 26-feature columns, with the same homogeneous scale.

## Retain a useful geometric statistic without fixing the new weights

Prokhorov's proved distance statistic is h(d)=1 for d<=7, h(d)=1/3
for 8<=d<=15, and zero otherwise, including disconnected nodes. Let
S be its sum over all (edge, containing three-face) flags. Then

    S=sum_{edges e} (2 f_7(E_e)+f_15(E_e))/3,
    S+K_3+B >= 7T.

There is an exact endpoint marginal, so S need not remain a free edge
port. Write W(E)=(2f_7(E)+f_15(E))/3. For a finite vertex forest F,
its contribution to S is

    S_v=(1/2) sum_{u in F} W(F minus u)
       =(1/2) sum_{pairs {a,b} in F} (n-d(a,b)-1) h(d(a,b)).

For a pair in a spherical tree, exactly n-d-1 node deletions leave its
unique path intact. This proves the second identity, with disconnected
pairs contributing zero. The factor 1/2 counts the two endpoints of
each actual edge. For an ideal vertex whose affine factors are T_i of
rank r_i, set H_i=sum_{u in T_i} W(T_i minus u). Its contribution is

    S_v=(V/2) sum_i H_i/(r_i+1).

A link vertex selects one node deletion in each affine factor; summing
over all selections gives this formula, including spectator factors.
Thus S=sum_v S_v is bound to the same vertex counts as P,Q,A and the
old features. This identity uses actual endpoint diagrams, not an
assumption that the endpoint diagrams agree beyond their shared edge.

For r>=17, the closed path and D formulas are

    W(A_r)=(29r-176)/3,
    J(A_r)=sum_{pairs}(d+1)h(d)=(205r-1696)/3,
    W(D_r)=W(A_r)+1,       J(D_r)=J(A_r)+3.

Each finite component contributes (n W-J)/2 to S_v. For 2D_17+D_23,
W=378, J=2208 and S_v=9669. These are Prokhorov's different distance
statistics, not the manuscript's ordinary-pair table sums 113.5 and
145. For every rank the closed histogram formulas are N_d(A_r)=r-d
for 1<=d<=r-1, and
N_d(D_r)=max(r-2-d,0)+2+[d=2] for 1<=d<=r-2; truncate the weighted
sum at d=15, not the allowed ranks. For an affine A factor with m
nodes, H=m W(A_{m-1}); the other affine families still need their
own symbolic coverage in the full tail implementation.

Here E_e is the spherical edge diagram of rank n-1 and f_l counts
unordered node pairs at distance at most l. These are fixed geometric
data and a source-backed necessary inequality. Their dual multipliers
are free. Adding this row does not fix the manuscript's corner-weight
coefficients. Replacing h itself by a richer weight law is another
possible extension, requiring shared three-face states as above.

Together with the two apex/link bounds, this implies

    S+K_3+A/2 >= 7T,       S+K_3+Q/3 >= 7T.

For example, the first inequality enters a finite vertex budget through
A_v/2. This is the controlled incidence compensation that cannot be
obtained by retuning the old independent budget alone.

## Tail validity, aggregation limits and a printed sign

The parity-specific Khovanskii comparisons on pp. 405–406 are

| Quantity | even n | odd n |
| --- | --- | --- |
| E_3/T | <12(n-1)/(n-4) | <12n/(n-3) |
| K_3/T | <6(n-2)/(n-4) | <6(n-1)/(n-3) |
| K_4/f_4 | <8(n-3)/(n-6) | <8(n-2)/(n-5) |

Use their weak closures in an outer relaxation. Clearing the positive
denominators for n>=8 gives affine polynomial coefficients in n.
All the new incidence and capacity rows hold uniformly; they require
no maximum diagram rank, dimension or facet count. Split the integer
tail into n=N_odd+2t and n=N_even+2t, t>=0, when forming tail SOS
identities. A finite reduction of the *whole* full-feature model still
requires coverage of its vertex, edge and shared-star parameters.

The last displayed odd-n line on p. 406 has a printed plus sign before
the bad-face fraction. The preceding inequality S/T>=7-K_3/T-B/T
requires a minus there; the following printed polynomial uses the
minus and is consistent. Derive the rows from the incidence argument.
That polynomial is n^2-(96C+68)n+147<0. At n=57 and C=0 its left side
is -480. Thus even the unsimplified odd-dimensional comparison cannot
exclude 57 by changing only the nonnegative edge constant C.

The stronger caps at k=6,7,9 do not by themselves improve the
bad-face fraction used in Lemma 1 near dimension 57. The line
b<=2k-10 is attained by the permitted (k,b)=(8,6),(10,10) states, so
their mixtures already attain its upper envelope for average k in
[8,10]. Those are relaxed states, not assertions of realization.
Link/Gram incidences must prevent these mixtures or control their
vertex costs before a stronger conclusion follows.

## Integration in the global formulation

The [incidence module](../../research/global_formulation/incidence.py)
retains these uniform formulas. The
[global builder](../../research/global_formulation/model.py) binds its
P,Q,A,S ports to the same whole-vertex measures that supply all 26
manuscript features. Shared face Gram schemas cover every triangle
and quadrilateral. Their marginal equations and the actual-polytope
lift are specified in [global_formulation.md](global_formulation.md).
The dual finds all corner coefficients together with the new multipliers.

The stronger Andreev incidence strata and overlaps between different
face stars remain possible refinements; they are not all imposed by
this initial sparse moment level. Source review and exact certificate
verification are still required before accepting an improved bound.
