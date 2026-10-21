# Large A components with arbitrary spectators: transfer of the whole-A certificate

20 October 2026. Complete paper candidate, awaiting independent review.
No scientific execution, optimizer, lattice classification, Lean or TeX.
Accepted sources are unchanged. No arithmeticity is assumed.

**Theorem.** A vertex farthest from a strict interior point of a compact
hyperbolic Coxeter polytope cannot have a WHOLE A_r COMPONENT with
21<=r<=29, regardless of its other finite components. Together with the
already adopted universal compact dimension ceiling29, this excludes
every relevant A component of rank at least21 at a farthest vertex.
It does not give a new unconditional ambient dimension ceiling.

The new points are the zero-own-contact classification, exact clearing
of spectator pressure by ACTUAL feasible points, and the fact that every
weak central opposite has no spectator contact. No metric face is
assumed Coxeter, and no projected wall is assumed actual or native.

## 1. Pinned inputs and notation

All source paths below are under `results/`.

* Whole-A21 proof:
  `compact_whole_a21_farthest_vertex_exclusion_candidate_20261020.md`, SHA256
  `3abe11bb70ce677069279462a91111ed6b8ab8889e3ec9e656ac71c1094200fb`.
* Whole-A21--29 extension:
  `compact_whole_an_geometric_a21_plateau_transfer_20261020.md`, SHA256
  `8a58b7747e1a3d8c9affe7df17038781a62983ddd70ccb670502c6f40b1f9cee`.
* Actual finite-face pressure:
  `compact_farthest_spherical_face_projection_pressure_20261020.md`, SHA256
  `abfcb96b6eb74aefdc6498cda362d81b1d4a1f9c32ad17fe3b464cff26fba565`.
* The mixed actual-point and orthogonal-opposite inputs respectively have
  SHA256 `3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`
  and `4216d06c2e03fd336898b8dac289b9969fa096c8a1d08bd67ddeb8ab5ec2f23d`.
  Their named sources are `compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`
  and `compact_farthest_vertex_orthogonal_component_opposite_obstruction_20261020.md`.

Put N=r+1, 22<=N<=30. The component unit normals have duals

    D_i=omega_i^2=2i(N-i)/N,
    C_ij=2i(N-j)/N (i<j),
    h_i=x_A.omega_i, h_0=h_N=0,
    ||x_A||^2=(1/2)sum(h_(i+1)-h_i)^2<1.

The heights h_i are strictly concave, since the comparison point is
strictly interior to every old facet. For central indices5<=i<=N-5,
call i weak when h_i<=q_i-.61, q_i=sqrt(D_i). The accepted range is

    D_i>=85/11, 5/2<q_i<31/8,
    F(q_i)/q_i<761/500,
    F(q)=[1+(1-.61/q)^2]/[1-(1-.61/q)^2].          (1)

For an actual opposite c write

    c=A p-v_A-v_S, A>0,
    Q=||v_S||^2, A^2=||v_A||^2+Q-1.

Here S is the entire orthogonal finite spectator forest. Its projection
v_S equals the old-normal combination w_S=B_S^-1 beta_S, with
nonnegative coefficients. The actual endpoint retains S+c, so Q<1.
Within each spectator component c has at most one contact, since two
create a nonpositive cycle. Consequently Q>0 implies Q>=1/4, and
each contacted component costs at least1/4. These facts do not restrict
spectator rank or type.

## 2. General isolated-cap pressure clearing ALL spectators

Suppose c has just one component contact, -c.n_i=gamma>0, so
v_A=gamma*omega_i. Write d0=1-Q>0 and u=d0/gamma^2. Then

    A^2=gamma^2(D_i-u).

If u>1, the ACTUAL finite set S+n_i+c has Schur residual
d0-gamma^2>0. Applying the pinned finite-face pressure and clearing
all spectators gives the PURE component inequality

    x_A.gamma*(omega_i-n_i)>=H(A,d0-gamma^2),
    H(A,delta)=A/[1+sqrt(delta/(A^2+delta))].

Since x_A.n_i>0, setting B=sqrt(D_i-1) and t=sqrt(u-1) gives

    h_i>sqrt(B^2-t^2)/(1+t/B)>B-t.                 (2)

The final inequality follows by positive squares: its squared margin
against (B-t)^2 is (B-t)*t^2/(B+t)>0.

If 0<u<1, use instead the explicit actual combination

    xi=c+w_S+(d0/gamma)n_i.

It has xi.S=0, xi.c=0, xi.n_i=-gamma+d0/gamma<0 and

    xi^2=-d0*(1-u)<0.

Its p coefficient A is positive. Every other actual wall pairs
nonpositively with each normal in this positive combination, so xi
is a future timelike feasible point of the FULL polytope. Farthest
comparison, with no spectator term left, gives

    h_i>sqrt(D_i-u)-sqrt(u(1-u)).                 (3)

The strictness uses x_A.n_i>0 in the unabridged pressure vector
`gamma*omega_i-(d0/gamma)n_i`. If u=1 the same combination is a
nonzero future null feasible ray, so compactness excludes equality.
For example normalized p+t*xi gives points escaping to that ideal ray.

These arguments do NOT regard c+w_S as an actual facet normal, nor
assume its normalized angle with n_i is a native Coxeter angle.

### Every isolated central opposite is strong

The local affine witnesses excluding isolated own labels3,4,6 and
parallel contact in the whole proof remain principal diagrams with
spectators present. Thus gamma is cos(pi/5), cos(pi/m) with m>=7,
or dotted gamma>1. If Q=0 all accepted isolated strong-floor proofs
apply verbatim inside the p+A span.

If Q>0, gamma>=phi/2, d0<=3/4, and therefore

    u<=3/phi^2=6-3phi,
    u-1<=5-3phi<73/500<(5/12)^2,

where phi=(1+sqrt5)/2>809/500. For u>1, (2), q_i>11/4 and
sqrt(D_i-1)>5/2 show

    q_i-h_i<4/21+5/12=17/28<.61.

For u<1, (3) gives

    q_i-h_i<(4/21)u+sqrt(u(1-u))
        <=2/21+sqrt457/42<127/210<.61.

The last radical uses sqrt457<107/5. Thus arbitrary foreign pressure
does not produce an isolated weak cap: every isolated central cap is
strong in the ACTUAL component projection.

## 3. Complete central old-component inventory and zero own contact

At the actual opposite of i, the retained old component forest is
A_(i-1)+A_(r-i), both paths of rank at least4. Restricting its actual
positive-definite replacement to these paths and c gives exactly the
accepted whole-A possibilities: no white contact; one white3 at an
admissible A/D/E position; one end white4; or two end white3 contacts,
one on each path. The two-path connected rank is r>=21, so its only
finite completion is A_r: an internal branch cannot have the two
short arms of D, a4 is not a path-end bond, and E/F/H exceptions are
too small. Extra spectator contacts cannot enlarge this list.

Unlike the whole-vertex proof, gamma=-c.n_i=0 cannot be dismissed by
assuming the entire old star is finite. Here is the complete repair.

No component contact at all has zero outgoing derivative and cannot
close i. Two white3 contacts give an exact all-3 affine cycle in the
full A_r+c diagram. For a single white3 at a global internal node
3<=j<=r-2, its three arms have lengths1,j-1,r-j. They contain affine
E8 arms1,2,5 or affine E7 arms1,3,3. A single white4 at a global
internal node gives the exact affine star3,3,4. Spectators do not
change these principal diagrams.

Only the global far endpoints1,r and next-to-end nodes2,r-1 remain.
Let m>=5 be the relevant retained path length plus one, measured from
that far endpoint; the other side is symmetric.

* A3 at node1: own full cost is (N-1)/(2N), retained white cost is
  (m-1)/(2m), and the height square is less than1/(2m)<=1/10.
  Its positivity requires a spectator contact. Any such component
  costs at least1/4, giving a reciprocal shadow above one in that
  orthogonal component, contradicting orthogonal-opposite obstruction.
* A3 at node2: own cost is (N-2)/N<1 and Q<2/m<=2/5. At most one
  spectator component can be contacted. If present its full cost Q
  exceeds A^2=Q-2/N, again forcing the same opposite to close a
  direction of an orthogonal component. If absent the height is
  negative.
* A4 at node1: retained white cost is (m-1)/m, so Q<1/m<=1/5;
  hence Q=0 and its own cost (N-1)/N gives negative height.

In the shadow arguments the contacted inverse diagonal gives
`T_j^2>=q_j/A^2>1`; all spectator components are retained in the
actual finite endpoint. This proves gamma>0 for every central cap
without a whole-star or harmless-spectator assumption.

## 4. Foreign contacts make every remaining exceptional role strong

The actual edge formula is

    cosh L_i=sigma*q_i/sqrt(delta0-Q),             (4)

where sigma and delta0 depend only on the old-component contacts.
Thus adding Q only strengthens any edge-only lower bound used in the
whole proof. All its nonexceptional roles remain strong. Its local
affine exclusions of own labels3,4,6 likewise remain unchanged.
The only roles requiring a new check are far-end3 and near4/simple3.

For a far-end3 helper, parameterize its white path by m>=5. Then
`delta0=(m+1)/(2m)` and `sigma=gamma+1/(2m)`, with
gamma>=cos(pi/5)>f=101/125. If Q>0, Q>=1/4 gives

    sigma/sqrt(delta0-Q)
       >(2fm+1)/sqrt(m(m+2))
       >=227/(25sqrt35)>761/500.                 (5)

The middle function is increasing for m>=5: its derivative has sign
`(2f-1)m-1>0`. The final squared numerator is
`4540^2-35*761^2=342365>0`.

For the near4/simple3 role, sigma>=1/sqrt2+(m-1)/(2m)>11/10 and
delta0-Q<=(m+2)/(4m)<=7/20. Thus its squared ratio is greater than
121/35>(761/500)^2. Equation(1) and (4) make both foreign roles
strong. The conclusion is therefore exact:

    EVERY weak central cap has Q=0.

Its complete profile and its height now coincide with those of the
accepted whole-A proof: a near4/simple3 packet, or the literal far
helper `v_A=gamma*omega_i+omega_e/2`, e=1 orr, with native own5,
native own m>=7 or a dotted own contact. No white or foreign demand
has been omitted to obtain this conclusion.

## 5. Transfer of the complete geometric certificate

All remaining steps of the pinned A21--29 proof now take place in the
same p+A span, with the same actual constraints:

1. A weak near4 cap has Q=0, so the accepted actual mate and its
   normalized height>sqrt3/2 are unchanged. The packet exclusion and
   opposite strong-anchor Cauchy argument remove weak near4 packets.
2. All other weak caps have Q=0. Their mutual products, heights and
   native gaps are exactly the accepted ones. The two-helper weighted
   AM--GM/native-product argument and finite mixed replacement thus
   exclude adjacent weak literal packets with the same constants.
3. Every selected finite or timelike pressure point in the weak5
   table is a positive acute construction from these actual caps and
   old A normals. Adding all old spectator normals is an orthogonal
   positive-definite block; those normals pair zero with the selected
   points. Their feasibility against ALL other actual walls follows
   from the original positive-combination argument. All pressure
   vectors, time coefficients and exact positive-dual certificates
   therefore remain unchanged.
4. Strong central floors, strict concavity and the norm bound all
   hold for x_A itself. The exact plateau lift, endpoint removal and
   six-hole energy estimate consequently give

       ||x_A||^2>4013/4000>1,

   contradicting ||x_A||<1.

For finite mixed replacements, appending the spectator block restores
the complete retained vertex rank without changing its Gram inverse
on the selected span. There is no inference that an arbitrary metric
face is Coxeter or has the full integral reflection arrangement.

This transfers the COMPLETE accepted certificate, not just its scalar
edge floors. It proves the claimed component exclusion for ranks21--29.
The global farthest-diagram alternatives of smaller A components and
other large components remain open; no new dimension ceiling follows
solely from this source.

## Contributions and review boundary

Finite57 derived the spectator-clearing actual points, zero-own-contact
repair and foreign-weak-cap exclusion and wrote this transfer. Root and
guard authored or co-developed several pinned whole-A and actual-point
inputs. Guard independently checked the generic clearing construction
in messages while investigating the separate B lane. This is a new
collaborative proof candidate, not an outside-author review.
