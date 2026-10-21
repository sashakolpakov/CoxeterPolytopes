# Sharp reciprocal transport and4462 ideal-six vertices

For an entirely ideal finite-volume right-angled hyperbolic six-polytope,

    V >= ceil(2673F/136), and therefore V >=4462 when F>=227.

Here V counts ideal vertices and F facets. This strengthens the original
4425 bound by a sharp estimate of the formal cusp-frame reciprocal average.
It leaves the dimension bounds and existence question unchanged. It uses
the universal cusp geometry and face counts in
[the original transport proof](ideal_ra6_cusp_transport_research_20261018.md),
whose separate [independent review](ideal_ra6_cusp_transport_independent_review_20261018.md)
passes its original4425 conclusion. Those original files are unchanged.

## Exact maximizer of the reciprocal average

Let epsilon_1,...,epsilon_5 be independent uniform signs, mu=3/2, and
sum p_i^2=1/4. Then the sharp inequality is

    E_epsilon [1/(mu-sum epsilon_i p_i)] <=103/132.                  (1)

Equality occurs precisely when p_i^2=1/20 for every i. This is a maximum
of formal frame averages. It is not an assertion that every formal
neighbor is an actual edge or that a polytope attains the resulting bound.

Signs of p_i can be absorbed in the independent epsilon_i. Put a_i=p_i^2,
so a_i>=0 and sum a_i=1/4. This is a compact simplex. For
X=sum sqrt(a_i)epsilon_i, Cauchy gives

    |X| <=sum sqrt(a_i) <=sqrt(5)/2 <mu.

Consequently the following positive even-moment series converges uniformly
over the whole simplex, and defines a continuous function there:

    E[1/(mu-X)] = (1/mu) sum_(k>=0) E[X^(2k)]/mu^(2k).             (2)

Symmetry kills the odd terms. Uniform convergence follows from the
geometric majorant (sqrt(5)/3)^j for the unsymmetrized series. In
particular the maximum exists, including at boundary points with zero
coefficients.

Fix two squares a,b and leave all other coefficients unchanged. Write
X=R+Y, where R is the independent symmetric sum of the other three terms
and Y=sqrt(a)epsilon+sqrt(b)eta. If s=a+b and t=ab, then

    E[Y^(2j)] = E[(s+2sqrt(t)epsilon eta)^j]
               =sum_(0<=l<=floor(j/2)) binom(j,2l) s^(j-2l)(4t)^l. (3)

At fixed s every coefficient is nonnegative, so replacing a,b by s/2,s/2
increases or preserves every even moment of Y. Expanding against R gives

    E[(R+Y)^(2k)]
      =sum_(j=0)^k binom(2k,2j) E[R^(2k-2j)] E[Y^(2j)],

again a sum of nonnegative terms. Thus balancing a,b cannot decrease any
even moment of X. It increases the fourth moment by exactly

    4(s^2/4-ab)=(a-b)^2.

If a!=b, the positive fourth-moment coefficient in (2) makes the reciprocal
average strictly larger after balancing. Therefore no global maximizer
can have unequal squares. Every maximum has a_i=1/20. This is a global
maximizer proof; no finite list of sampled coefficient directions is used.

At equal squares, X=(sum epsilon_i)/(2sqrt(5)). The squared values5/4,
9/20 and1/20 have multiplicities2,10 and20 among the32 sign vectors.
Pairing opposite sign vectors gives exactly

    E[1/(mu-X)] =mu E[1/(mu^2-X^2)]
      =(3/2)[(1/16)/1 +(5/16)/(9/5) +(5/8)/(11/5)]
      =103/132.                                                    (4)

Equations(2)--(4) prove (1), including equality and the boundary of the
coefficient simplex.

## Consequences for actual cusps and faces

The reviewed transport proof canonically normalizes every cubic cusp and
preserves h=z-v/8 across good edges. Choosing future null Z perpendicular
to h gives E(v)=-Z.v>0 and p_i=(Z.e_i)/E(v), with sum p_i^2=1/4. The
neighbor ratio is r=3/2-sum epsilon_i p_i. All32 formal ratios are positive,
so adding nonexistent neighbors is an upper relaxation. Equation(1) gives
a ratio sum at most32*103/132=824/33. With f=1/E, each actual undirected
good edge contributes a reciprocal pair f(u)/f(v)+f(v)/f(u)>=2. The
unweighted total good degree on cubic cusps is therefore at most
(824/33)V_c.

At each six-cusp the number t of triangular link edges satisfies
t<=64+g/2. Hence the cubic total is at most

    [64+(824/33)/2]V_c =(2524/33)V_c.

The reviewed noncubic pointwise bound72 is smaller than2524/33. Since
each triangular two-face has three ideal vertices, the total number T
of triangular two-faces satisfies

    3T <= (2524/33)V.

The reviewed ideal-six identities are a_2=25V+a_4-F and
T>=4a_2-80V=20V+4(a_4-F). Multiplying the latter by three before combining
with the corner bound yields

    60V+12(a_4-F) <=(2524/33)V,
    a_4-F <= (136/99)V.

Every ideal five-facet has at least56 four-facets by the earlier universal
carrier theorem, so a_4>=28F. Therefore

    27F <=(136/99)V,
    V >=(2673/136)F.

The inherited F>=227 gives2673*227=606771 and
V>=ceil(606771/136)=4462. No triangular corner, good-edge endpoint or
facet incidence is counted with a new convention here.

## Exact arithmetic and scope

[The tiny producer](../scripts/fresh_ideal_ra6_sharp_reciprocal_20261018.py)
and [certificate](ideal_ra6_sharp_reciprocal_20261018.json) check the equal-sign
sum and every rational consequence above. They also compare coefficients
in (3) with the direct pair expansion through j=16. Those finite coefficient
checks are checks of the formula, not the universal all-moment proof;
the latter is the identity Y^2=s+2sqrt(t)epsilon eta and the written
compact-maximizer argument. The source imports no cusp-model producer.

Generation and fresh saved replay exited0/PASS in chunks f005b4 and76467e,
with peak RSS26,771,456 and26,525,696 bytes respectively. Both held the shared
exclusive lock and used the59s alarm, start/end resource gates and1GiB RSS
check. No process remains. No Lean, large result
load, optimization, new example or dimension exclusion is claimed.
The original4425 source/certificate and all earlier frozen evidence
remain byte-identical.
