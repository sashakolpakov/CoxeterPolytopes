# Two actual logarithmic transports and finite-edge pairing in dimension ten

20 October 2026. Paper inequalities and a prepared relaxation, not a new
dimension theorem. No scientific job, optimizer, census, checker, Lean or
build was run in preparing this note. The two logarithmic estimates were
proposed independently by the root coordinator and checked here.

The geometric source for the second move is
`results/right_angled_r7_four_direction_transport_research_20261020.md`
(SHA256 `0407d5ddc6d0525fb18be33e86821a83cebfdc42ce8479aae495d20f5fc362b4`).
The actual P3 construction is in `manuscripts/right_angled_dimension10.tex`.
The upper defect input is
`results/right_angled_three_face_upper_d_two_excess_20261020.md`
(SHA256 `b7b8d841c03fd813d31a53853255429803fdc6ec44c6103703e138a708423f8f`).
These are explicit written geometric inputs; no finite computation below
substitutes for their actual-wall and reversal proofs.

## 1. One null vector per actual component

At a finite vertex q of a hypothetical finite-volume right-angled
ten-polytope, let n1,...,n10 be its outward orthonormal normals. A direction
is ideal-ending if its actual edge ends at a cusp. Write k for their number
and m=10-k for the number of finite-ending directions.

The P3 graph has the actual opposite

    qS=2q-sum_(i in S) ni, |S|=3,
    h=3q-sum_i ni, h^2=1.

It is present precisely for a P3 three-face. On each actual finite graph
component choose one future null Z perpendicular to its invariant h. Put
F=-Z.q and p_i=-Z.n_i/F. Then

    sum p_i=3, sum p_i^2=1, rS=F(qS)/F(q)=2-sum_S p_i.

The R7 graph is a different actual graph. A four-set A of finite-ending
directions gives a move if each of its four triples is central R7. Its
actual opposite and invariant are

    qA=3q-sqrt2 sum_A ni,
    h4=2sqrt2 q-sum_i ni, h4^2=2.

The complete-wall feasibility and reverse R7 identification are proved in
the pinned source above. On each of these components choose its own single
future null Z perpendicular to h4. Its coordinates satisfy

    sum p_i=2sqrt2, sum p_i^2=1,
    rA=F(qA)/F(q)=3-sqrt2 sum_A p_i.

No identification of the two null vectors or their coordinates is used.
Each graph separately is finite and undirected. On each undirected edge
the two directed logarithms are opposite. Therefore in BOTH graphs

    sum_(actual directed moves) -log(F'/F)=0.                 (1)

This is an exact ordinary vertex sum, not a weighted-height average.

## 2. A reverse Jensen bound with its scalar proof

If a positive random variable r is bounded below by a>0, has mean mu and
variance V, then

    E[-log r] <= -log mu + V/(2*a*mu).                       (2)

Indeed for x=r/mu<=1 the scalar remainder estimate with denominator
2*r*mu is equivalent to

    1-x^2+2x log x >=0.

Its derivative is 2(1-x+log x)<=0 and its value at one is zero. For x>=1,
`x-1-log x <= (x-1)^2/2`, as follows by differentiation. Since a<=r and
a<=mu, the same upper coefficient 1/(2*a*mu) works on both sides.
Take the expectation; the first-order term vanishes.

Suppose now mu=mu0+d*delta, delta>=0, and

    V<=A*delta-B*delta^2,  B>=0,  A/d<=2a.

Then (2), discarding its negative quadratic variance term, and
`log(mu/mu0)>=d*delta/mu` give

    E[-log r] <= -log mu0.                                  (3)

This remains valid when delta=0. The sharp constant is attained by the
uniform selected coordinates and uniform complementary coordinates below.

## 3. Exact selected-subset variance

For ten real coordinates with sum b and squared norm one, select a fixed
j-element set I, j<10, and write S=sum_I p_i. Cauchy gives

    S<=M=(j*b+sqrt(j*(10-j)*(10-b^2)))/10.

Put delta=M-S and W=sum_I(p_i-S/j)^2. The complementary Cauchy bound is

    W<=1-S^2/j-(b-S)^2/(10-j)
      =2sqrt(10-b^2)*delta/sqrt(j*(10-j))
       -10*delta^2/(j*(10-j)).                              (4)

For a uniformly chosen l-subset T of I,

    Var(sum_T p_i)=l*(j-l)*W/(j*(j-1)).                     (5)

Equation (5) follows by expanding the diagonal and off-diagonal inclusion
probabilities l/j and l(l-1)/(j(j-1)); it requires no independence.
Equations (4),(5) apply to arbitrary real coordinates, including negative
ones, and to a fixed actual direction set I. They do not replace actual
move families by all subsets: unavailable moves are removed explicitly.

### P3: sharp selected sizes 3 through 9

Take b=3,l=3 and ratio r=2-sum_T p_i. Its full formal range is

    (11-sqrt21)/10 <= r <= (11+sqrt21)/10,
    16/25 < r < 25/16.

For an ideal-ending direction set of size k=3,...,9,

    mu0(k)=11/10-(3/10)*sqrt((10-k)/k), d=3/k,
    A/d=2*(k-3)/((k-1)*sqrt(k*(10-k))) <=1/2<2*(16/25).

The finite inequality, after squaring positive sides, has differences

    (k-1)^2*k*(10-k)-16*(k-3)^2
      =84,200,336,456,500,384,0   (k=3,...,9).

Thus define

    Lk=-binom(k,3)*log(mu0(k)), 3<=k<=9;  L0=L1=L2=0.        (6)

These are sharp upper bounds for the complete formal selected sums.
In particular L9=0. For k=10 use the full mean11/10 and variance7/300
in (2), with a=16/25 and log(11/10)>=2/21:

    L10 <=120*(-2/21+35/2112)=-5815/616.                    (7)

### R7: sharp selected sizes 4 through 9

Take b=2sqrt2,l=4 and ratio r=3-sqrt2 sum_T p_i. Its range is

    (7-2sqrt6)/5 <=r<=(7+2sqrt6)/5,
    21/50<r<12/5.

For a finite-ending direction set of size m=4,...,9,

    mu0(m)=7/5-(4/5)*sqrt((10-m)/m), d=4sqrt2/m,
    A/d=4*(m-4)/((m-1)*sqrt(m*(10-m))) <=5/6<2*(21/50).

The positive-side squared differences are

    25*(m-1)^2*m*(10-m)-576*(m-4)^2
      =5400,9424,12696,13716,10384,0   (m=4,...,9).

Set

    Mm=-binom(m,4)*log(mu0(m)), 4<=m<=9;  M0=...=M3=0.       (8)

In particular M8=0 and M9=-126log(17/15). For m=10, the complete mean
is7/5 and variance8/75. Equation (2), a=21/50 and log(7/5)>=1/3 give

    M10<=210*(-1/3+40/441)=-1070/21.                        (9)

Both logarithm bounds use `log(1+x)>=2x/(2+x)` for x>=0; its derivative
comparison has nonnegative numerator x^2. These strengthen the first
paper estimates -785/88 and -860/21 without altering any geometric input.

## 4. Ordinary errors, defect and the two global rows

Let f_k count finite vertices of ideal degree k, and f=sum f_k. Let B be
the number of P3 three-faces. Each P3 has two finite corners, so

    H=sum_k binom(k,3) f_k-2B >=0.

It counts exactly the missing P3 moves among ideal-ending triples. The
reviewed three-face argument gives H<=2E, with

    E=sum_(Y!=P3)(facets(Y)+cusps(Y)-10)>=0.

Let z_k count noncentral finite-direction triples at ideal-degree-k
vertices; these are exactly the big-three-face A0 incidences. Put Z=sum z_k.
For m=10-k, the number b4 of unavailable four-moves obeys

    b4<=sum_k (m-3)_+ z_k<=7Z,
    0<=z_k<=binom(m,3)f_k,
    4Z<=10E-D,   -30E<=D<=2E.                            (10)

Here D=sum(A2-3A0) is an ordinary count over all actual three-faces.
The small types have defect zero. Its equivalent vertex identity is

    D=sum_k [binom(k,2)*(10-k)-3binom(10-k,3)] f_k.          (11)

Removing a formal P3 logarithm costs at most log(25/16); removing a
formal R7 logarithm costs at most log(12/5). Equation (1) therefore gives
the simultaneous necessary inequalities

    0<=sum_k Lk f_k+log(25/16)*H,                          (12)
    0<=sum_k M_(10-k) f_k
       +log(12/5)*sum_k (7-k)_+ z_k.                      (13)

The weaker replacement of the last error by `(7/4)log(12/5)(10E-D)`
is safe. The two rows do not require a common null vector between graphs.

## 5. Actual finite-edge pairing

Call a finite--finite edge regular if every three-face containing it is
one of P3,Q8,R7. Such an edge is never contained in a P3. The Q8 and R7
incidences force every transverse status to reverse: a pair of finite
directions is central R7, a pair of ideal directions is a leaf of R7,
and a mixed pair belongs to Q8. Since one of the two transverse classes
has at least two members, these pair statements fix each individual
status, including when the smaller class has size zero or one. Therefore

    k(q')=9-k(q).

Let N_k count irregular finite-edge incidences at degree-k vertices. Exact
pairing and the big-three-face union bound give

    (10-k)f_k-N_k=(k+1)f_(9-k)-N_(9-k), 0<=k<=9,
    0<=N_k<=(10-k)f_k,
    sum N_k<=30E,
    sum N_k<=20E-D-2Z.                                   (14)

The last bound uses `2e_FF=2p-d-2A0-2A3` on each big three-face and
p<=10e. There is no extra factor two in the sum of absolute discrepancies
over the five unordered pairs. These are actual edge counts, independent
of both diagonal opposition graphs.

## 6. Dimension-ten ordinary identities and the prepared relaxation

Let c count cusps and a_j count j-faces. The cusp three-face incidence is
2^7*binom(9,2)=4608. Let R=sum_(Y!=P3)(8-facets(Y)); **R is signed**.
The ordinary identities are

    a3=4a2-30f-2304c,
    E=B-10a3+8a2+4608c,
    R=8a3-8a2-2B
     =-5B/4-15f+2304c-3E/4.                              (15)

The reviewed exceptional-corner bound remains X<=18E+4R and

    X>=sum_k [binom(10-k,3)+binom(k,2)*(10-k)] f_k.         (16)

It is incorrect to replace R by a nonnegative polygon excess. Likewise,
the dimension-eleven weighted row cannot simply be scaled down to this
dimension. A preparatory LP may include the actual dimension-ten
truncation Dehn--Sommerville equations and inherited facet minima, but a
floating-point infeasibility status would only nominate an exact dual.
Conversely a feasible relaxation is not a polytope construction.

The accompanying scout source is prepared separately for root review and
serialized execution only. No new dimensional conclusion is made here.

The separately reviewed ideal dimension-four theorem supplies additional
finite-face cover rows, proved with a finite-vertex parity seed in
`results/right_angled_finite_vertex_seed_and_face_cover_20261020.md`.
These are not consequences of the logarithmic rows alone. The prepared
relaxation explicitly includes and source-binds this new geometric input.
