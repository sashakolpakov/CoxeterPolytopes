# A common-white weighted drift for compact Coxeter polytopes

20 October 2026. Paper proof only; no computation, dimension exclusion, or
new classification is asserted. This refines the undirected finite-edge
drift discussed in the current uniform-bound research. All weights below
are shared across an actual edge, rather than selected separately at its
two endpoints.

## 1. A global identity

Let P be a compact hyperbolic Coxeter n-polytope, with future unit vertex
vectors p_v and outward unit facet normals. Compact Coxeter polytopes are
simple, so the n normals at a vertex v have a positive definite Gram matrix
B_v. Write omega_(v,i) for their dual vectors and

    d_(v,i) = (B_v^(-1))_(ii) = ||omega_(v,i)||^2.

Choose ONE positive function f on the isomorphism classes of finite,
unit-normal Gram diagrams of rank n-1. An edge e=v--w has length L_e>0 and
common white diagram W_e. Set

    k_e = f(W_e) sqrt(det B_(W_e)) / sinh L_e.

This definition is undirected. At v put

    rho_f(v) = sum_i f(W_(v,i)) omega_(v,i),
    A_f(v)   = sum_i f(W_(v,i)) sqrt(d_(v,i)) tanh(L_(v,i)/2).

Then the exact vector identity is

    sum_(w adjacent v) k_(vw) (p_w-p_v)
      = sqrt(det B_v) [A_f(v) p_v-rho_f(v)].                 (1)

Indeed the inward unit tangent of the edge released at i is
-omega_(v,i)/sqrt(d_(v,i)), whence

    p_w = cosh L_e p_v-sinh L_e omega_(v,i)/sqrt(d_(v,i)).

The principal-cofactor identity gives
sqrt(det B_(W_e))=sqrt(det B_v) sqrt(d_(v,i)), and
(cosh L-1)/sinh L=tanh(L/2). Substitution proves (1).

Summing (1) over the finite vertex set gives zero on the left. Since
rho_f(v) is spacelike and perpendicular to p_v, the condition
A_f(v)>=||rho_f(v)|| makes every nonzero right-hand summand future causal.
Positive f makes rho_f(v) nonzero. A finite sum of nonzero future causal
vectors cannot vanish. Consequently SOME actual vertex satisfies

    A_f(v) < ||rho_f(v)||.                                (2)

This uses compactness, the actual edge graph, and actual edge lengths.
It does not follow from a formal collection of endpoint diagrams. No
arithmetic assumption is needed.

## 2. Exact optimization at a whole A_n vertex

Suppose a vertex diagram is the connected A_n diagram. Put N=n+1 and

    d_i = 2i(N-i)/N,        d_0=d_N=0,
    K_ij = (B^(-1))_ij/sqrt(d_i d_j).

For any local positive weights f_i define w_i=f_i sqrt(d_i). The ratio
relevant to a uniform lower bound on tanh(L_i/2) is

    ||sum_i f_i omega_i|| / sum_i f_i sqrt(d_i)
      = sqrt(w^T K w)/(sum_i w_i).

Its minimum over positive w is 1/sqrt(T_n), where

    f_i^* = sqrt(d_i)-(sqrt(d_(i-1))+sqrt(d_(i+1)))/2,
    w_i^* = sqrt(d_i) f_i^*,
    T_n   = sum_i w_i^*
          = sum_i d_i-sum_(i=1)^(n-1) sqrt(d_i d_(i+1)).   (3)

The square-root sequence is strictly concave, so all f_i^*>0.
Furthermore K^(-1)=diag(sqrt d) B diag(sqrt d), hence
w^*=K^(-1)1. Cauchy--Schwarz proves the stated minimum. One can also see
the equality directly:

    rho_(f^*) = sum_i sqrt(d_i) n_i,
    ||rho_(f^*)||^2 = T_n = sum_i w_i^*.

An equivalent exact expression is

    T_n = 1+sum_(i=1)^(n-1) (1-r_i)/(1+r_i),
    r_i = sqrt[i(N-i-1)/((i+1)(N-i))].                    (4)

To verify (4), K_ij for i<j is the product r_i...r_(j-1).
The inverse quadratic form of this correlation matrix is
x_1^2+sum_i (x_(i+1)-r_i x_i)^2/(1-r_i^2); evaluate at x=1.
Equivalently each summand is tanh((u_(i+1)-u_i)/4), where
u_i=log(i/(N-i)). Thus T_n=(1/2)log n+O(1). This improves the constant
asymptotic ratio 4/(pi sqrt 3) obtained from f_i=1.

## 3. How these weights can, and cannot, be used globally

At an A_n vertex the white diagram released at i is
A_(i-1)+A_(n-i), with rank-zero components omitted. Its diagram
isomorphism class determines i up to reversal. Formula (3) has exactly
that symmetry. It therefore defines a legitimate single function f on
these white classes, which can be extended positively to every other
white class before applying (1)--(2).

This does NOT make each other vertex locally optimal, and (2) does NOT
say that the resulting vertex is of type A_n. Every occurrence of the
same white diagram must retain its assigned value. For example the
connected white A_(n-1) also occurs at a B_n end and a D_n whisker; other
collisions occur at disconnected finite vertex diagrams. A dimension
argument must control all such vertices and all simultaneous consumers.
One may not optimize independently at the two ends of an edge.

A simple globally defined alternative, without an A-type case choice,
is

    f(W) = (det B_W)^(-3/2).

After cancelling a common vertex factor, its coefficients at EVERY
vertex are f_i=d_i^(-3/2), with relative edge weights w_i=1/d_i.
For A_n their sum is exactly the harmonic number H_n, since

    1/d_i = (1/2)(1/i+1/(N-i)).

## 4. The still-open quantitative obligation

For ANY chosen shared f, abbreviate w_i=f(W_i)sqrt(d_i),
W=sum_i w_i, and r=||rho_f||/W at the vertex supplied by (2).
If tau_0=tanh(L_0/2)>r, the total weight of directions with L_i<L_0 obeys

    sum_(L_i<L_0) w_i > W-||rho_f||/tau_0.

This follows by retaining only the long-edge terms in A_f. A single
opposite wall can close several directions, so a bound on the number of
distinct opposite walls is not automatically a bound on this weight.

The identity and optimization above are reusable global lemmas. A
material dimension exclusion still requires a universal bound on the
compatible short-edge weight, covering multi-contact and shared closing
walls, arithmetic root-length classes when used, and all other vertex
types sharing the selected white weights. No such bound is proved here.
