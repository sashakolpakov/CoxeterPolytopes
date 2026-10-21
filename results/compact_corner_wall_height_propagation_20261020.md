# A wall-height floor for all other closing walls

20 October 2026. Own paper lemma and one mixed-corner consequence.
No computation, dimension ceiling, TeX or Lean is asserted.

## 1. Height propagation

At an actual finite vertex p of a compact Coxeter polytope use outward
unit normals, finite Gram B, dual vectors omega_i, and coordinates

    c=a p-sum_i beta_i omega_i,     a>0, beta_i>=0

for an actual wall c not through p. Put t=B^(-1)beta and let J be the
union of the connected old components touched by beta. The inverse of
each connected finite Coxeter Gram is entrywise positive. Thus

    mu_c=min{1/2, (1/2)min_(j in J) t_j}>0.

EVERY distinct actual wall d not through p has height

    b=-p.d >= mu_c/a.                                  (1)

To prove this, write d=b p-sum_j gamma_j omega_j. If d contacts J, a
nonzero native contact gamma_j is at least1/2, so their spatial product
S=beta^T B^(-1)gamma is at least mu_c. Since c.d=S-ab<=0, (1) follows.
If d does not contact J, S=0. Distinct actual Coxeter walls with negative
product have product at most-1/2, hence ab>=1/2>=mu_c.

The restriction that d is NOT through p is essential. Incident old
walls have zero height and do not have nonnegative contact vectors.

If d closes a finite edge of length L, its common-white Schur norm
delta satisfies0<delta<=1 and b=sqrt(delta)sinh L. Therefore every
edge NOT closed by c has

    L>=asinh(mu_c/a).                                  (2)

In particular all edges in old components untouched by c obey (2),
since their rays have identically negative product with c.

This lemma does not assert a uniform positive mu_c over all finite
diagrams. Along a long A-chain it can decay like1/n, so the previously
constructed arbitrary-rank short-edge star remains a valid barrier.

## 2. A mixed A1/A2 profile is impossible at a high-dimensional maximum

Suppose the whole vertex diagram consists only of A1 and A2 components.
Assume an actual outside wall c has just two nonzero contacts: label3
to one A1 normal and label5 to one node of an A2 block. Set
phi=(1+sqrt5)/2. Its squared height and propagation floor are

    h=a^2=phi^2/3+1/4-1=(4phi-5)/12,
    mu_c=min{1/2,1/4,phi/6}=1/4.

The wall c closes three actual directions: the A1 edge and both A2
edges. For each replacement the n-wall Gram is positive definite
(the relevant rank-three diagrams are H3 or A2+A1). The acute projection
lemma therefore supplies the actual finite endpoint, strictly inside
every omitted wall. Another facet cannot truncate it first.

The A1 edge and the uncontacted A2-mate edge have respectively

    cosh L_(A1)=sqrt3*phi/2,
    cosh L_(mate)=phi^2/sqrt3.

Consequently their half-tangent SQUARES exceed1/6 and1/5, respectively.
The first comparison reduces to phi>121/75, and the second to
phi>19/12; both follow by evaluating x^2-x-1 at the stated rationals.

Every other closing wall has

    b^2>=1/(16h)>40/81.

For the strict rational inequality use phi<13/8 in
3/[4(4phi-5)]>40/81. Thus every other edge has

    tanh^2(L/2)>1/10,                                  (3)

because sinh^2 L=4tau^2/(1-tau^2)^2 takes the value40/81
at tau^2=1/10. Here “other edge” means an edge not closed by c;
in particular it includes every edge in an untouched old component.

Now choose p as the farthest vertex used in
`results/compact_farthest_vertex_energy_obstruction_20261020.md`.
The touched A1 component contributes energy greater than1/6; the
touched A2 block contributes at least the squared mate half-tangent,
hence greater than1/5. All remaining components are untouched by c.
By (3), each remaining A1 contributes more than1/10, and each A2 more
than(4/3)(1/10). Thus each remaining rank contributes more than1/15.
The total covering energy therefore exceeds

    1/6+1/5+(n-3)/15 = 11/30+(n-3)/15.

This is greater than1 for n>=13, contradicting the farthest-vertex
energy inequality. Hence this SPECIFIC actual mixed profile cannot
occur at a farthest vertex of an A1/A2 forest in dimension at least13.
The genuine three-dimensional [3,5,3] simplex realizes this local
profile; the conclusion uses the additional spectator directions and
their compulsory height floor. Other mixed profiles remain open.
