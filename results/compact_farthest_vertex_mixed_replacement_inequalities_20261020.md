# Actual mixed replacements at a farthest vertex

20 October 2026. Paper proof only. No computation, enumeration, build,
optimizer or formal verification was run. The mechanism was developed
in discussion with root. This note gives necessary inequalities, not a
dimension ceiling or a complete native-cap classification.

## 1. Positive acute projection gives an actual feasible point

Let P be a compact hyperbolic Coxeter polytope, with outward unit facet
normals. Distinct actual normals have nonpositive products. Let p be an
actual vertex and let m_1,...,m_k be linearly independent ACTUAL facet
normals whose Gram matrix G is positive definite. Put

    h_a = -p.m_a >= 0,    lambda = G^(-1) h,
    D = h^T G^(-1) h,
    q_0 = p + sum_a lambda_a m_a.                      (1)

Then lambda is entrywise nonnegative. One elementary proof is to write
lambda=lambda_+-lambda_- with disjoint nonnegative parts. If lambda_-
is nonzero, then

    lambda_-^T h
      =lambda_-^T G lambda_+ - lambda_-^T G lambda_- < 0,

because the off-diagonal entries of G are nonpositive. This contradicts
h>=0. This is the only inverse-positivity input needed below.

The selected-wall products of q_0 are zero. For every UNSELECTED actual
facet normal m,

    q_0.m = p.m + sum_a lambda_a (m_a.m) <= p.m <= 0.

Moreover q_0^2=-1-D and p.q_0=-1-D, so q_0 is future timelike. Thus

    q=q_0/sqrt(1+D)                                  (2)

is an ACTUAL point of P on all the selected walls. If k=n=dim(P), it
is an actual finite vertex. No assertion that an abstract Gram diagram
is realized has been made: the normals in (1) must already be facets
of P. Nor is any chosen facet or lower-dimensional face assumed to be
Coxeter intrinsically.

## 2. The additional farthest-point inequality

Fix a strict interior future unit point z and let p maximize -z.p over
the actual vertices. Write z=E(p-x), with x in p-perp and ||x||<1.
The same maximum bounds EVERY point q of P. Indeed a normalized positive
convex combination of future unit vertices represents q; the norm of
their unnormalized convex combination is at least one, so -z.q<=E.
Equivalently, hyperbolic metric balls are geodesically convex.

Write the actual point (2) in the form

    q=(p-v)/sqrt(1-R),    v in p-perp,
    R=||v||^2=D/(1+D)<1.

If D>0, comparison of distances from z gives

    x.v >= F(R):=1-sqrt(1-R).                         (3)

In particular

    ||x||^2 >= F(R)^2/R
             =(1-sqrt(1-R))/(1+sqrt(1-R)).             (4)

For exact algebraic comparisons, put s=x.v. Since s<1, (3) is equivalent
to

    R <= 2s-s^2.                                     (5)

When D=0 the replacement is p itself and (3) reads 0>=0.

These inequalities apply to mixed vertices involving several outside
walls, as well as to actual projected face points. They supplement the
incident-edge conditions x.u_i>=tanh(L_i/2). In particular those edge
conditions alone do not establish that p is farthest from z.

## 3. Isolated-cap specialization

At p let the old unit Gram be B, C=B^(-1), d_i=C_ii, and
u_i=omega_i/sqrt(d_i). Let J be a set of old directions with distinct
actual opposite walls c_i, each contacting ONLY old wall i. Write

    c_i=a_i p-b_i u_i,    b_i^2-a_i^2=1,
    r_i=a_i/b_i=tanh L_i in (0,1).

Let K_J=(u_i.u_j)_(i,j in J), which is positive definite. The cap Gram
is

    G_cap=diag(b_i) (K_J-r r^T) diag(b_i).

It is positive definite exactly when

    R_J=r^T K_J^(-1) r < 1.                           (6)

If (6) holds, combine these actual caps with the old walls indexed
outside J. The caps are orthogonal to all those retained old walls,
so the complete replacement Gram is positive definite. Section1 gives
an ACTUAL finite vertex q_J. Define

    v_J=sum_(i in J) (K_J^(-1)r)_i u_i.

The equations c_i.q_J=0 give u_i.v_J=r_i, and hence

    q_J=(p-v_J)/sqrt(1-R_J),
    x.v_J >= 1-sqrt(1-R_J).                           (7)

All positive definite cap subsets provide such inequalities, not just
singletons. For J={i}, (7) is precisely the usual incident-edge bound.
For larger J it is an additional actual-vertex constraint. No sign
assumption on the coefficients of K_J^(-1)r is required for (7).

For comparison, the separate acuteness/native tests remain

    c_i.c_j=(K_ij-r_i r_j)/sqrt((1-r_i^2)(1-r_j^2)),

which must be zero or a Coxeter nonright product at most -1/2.
These alphabet tests are NOT replaced by (7).

## 4. Combining several actual replacements

For any finite collection of valid replacements and nonnegative weights
w_J, put V=sum_J w_J v_J. Equation(7) implies

    x.V >= sum_J w_J F(R_J).

If some weight with R_J>0 is positive, then V is nonzero and

    [sum_J w_J F(R_J)]^2 < ||V||^2.                   (8)

The strict inequality uses ||x||<1 and Cauchy--Schwarz. A contradiction
to (8), established from complete actual native-cap data, would therefore
exclude the proposed farthest vertex. Establishing such a complete
cap/vertex coverage remains a separate obligation.

## 5. The previous B20 edge-only witness fails a two-cap replacement

The earlier formal whole-B20 relaxation used
tau_1=14/25, tau_2=63/100 and r_i=2tau_i/(1+tau_i^2), so

    k=K_12=1/sqrt(2),
    r_1=700/821,       r_2=12600/13969.

Its claimed vector x had x.u_i=tau_i, and squared norm less than
985816/10^6. The first coordinate bound there is non-strict, since
x_1=560/1000 exactly; the other bounds make the total strict.

For J={1,2}, compute

    R=2r_1^2+2r_2^2-2sqrt(2)r_1r_2,
    s=2(tau_1r_1+tau_2r_2)
       -sqrt(2)(tau_1r_2+tau_2r_1).

Direct rational bounds .8526<r_1<.8527, .9019<r_2<.9021 and
1.4142<sqrt(2)<1.4143 give

    .9048<R<.9069,      0<s<.618,
    2s-s^2 < .854076 < R.

Thus (5) fails with a strict margin greater than .05. This is an
explicit demonstration that the new actual-replacement requirement
goes beyond the previous pair/edge/energy relaxation. Those formal
caps did not satisfy the native Coxeter alphabet either; this subsection
does not present them as a Coxeter polytope or an actual counterexample.
It is a failed necessary-condition witness, not a new dimension bound.
