# Large B components: whole-B23 exclusion with arbitrary spectators

20 October 2026. Complete paper candidate, awaiting independent review.
No scientific execution, optimizer, Lean, TeX or accepted-source edit.
The foreign-far estimate is root's; the general isolated-clearing
construction was supplied by finite57 and independently checked here.

**Theorem.** At a farthest vertex of a compact Coxeter polytope, there
is no WHOLE B_r component with23<=r<=29, even with arbitrary orthogonal
finite spectator components. Combined with the already available
ambient compact ceiling29 this excludes all B_r components of rank>=23.
It is not itself a new ambient dimension bound or a rank/14 rate.

Use the accepted whole-B normalization:

    n_i=(e_i-e_(i+1))/sqrt(2) (i<r), n_r=e_r,
    omega_i=sqrt(2)*prefix_i, D_i=2i (i<r),
    h_i=x_B.omega_i, x1>...>xr>0, ||x_B||<1.

Call i weak if h_i<=sqrt(2i)-.61. The accepted whole-B sources are
`compact_whole_b_bulk_weak_packet_and_endpoint_pressure_20261020.md`,
SHA256 `e5a3c51858bb7b1e392b2828c874387be8d2044d5ae1be583a546b67eaddb54c`,
and `compact_whole_b_short_tail_extension_20261020.md`, SHA256
`b6adc03a5950a60a1e263aeb79869785268f1bb7905e74ec47a784446c989c5c`.
The proof below repairs precisely the inventory/point interfaces that
their WHOLE-vertex hypotheses did not cover.

## 1. Full finite white inventory and zero-own-contact repair

For5<=i<=r-4 the old white B part is A_(i-1)+B_(r-i), both ranks>=4.
Restricting the ACTUAL finite endpoint to these pieces and c gives
exactly the accepted B inventory: none; one left path contact
(end3/4, next3, or admissible E3); one right3 at LONG endpoint i+1;
or simple3 contacts joining the two pieces into B_r. Spectators can
add contacts, but cannot alter this finite principal restriction.
If both pieces are joined, c is an internal node of the connected
white B_r. A foreign contact would create a fork in a finite component
containing a4, impossible. Thus this branch has NO foreign contact.

All accepted affine principal witnesses remain exact independently
of contacts outside their chosen nodes. The sole old zero-own argument
that used positive index needs replacement: c may have simple3 at
far left1, zero at i, and foreign contacts. Write Q for their full
inverse-Gram cost. Then

    A^2=Q-1/2>0,
    Q<(i+1)/(2i), hence A^2<1/(2i)<=1/10.

The second line is the actual deletion-i finite residual. Any contacted
foreign node j has squared derivative
(C beta)_j^2/C_jj>=beta_j^2 C_jj>=1/4>A^2. The actual reciprocal-
shadow criterion therefore makes c close j too. This contradicts
the farthest orthogonal-component opposite theorem. All other
zero-own cases have the accepted affine fork, double4 path or cycle;
with no B contact the B derivative is zero. No foreign contact has
been silently removed.

For nonisolated profiles the accepted edge Schur numerator is
unchanged and the residual is reduced by Q. Thus every profile proved
strong by a direct edge estimate there remains strong. The only
unresolved nonisolated candidate is far3@1 plus own gamma>=cos(pi/5).

## 2. Every foreign isolated cap is strong

Suppose c has only its own contact gamma at B node i and foreign
spatial vector vS, with Q=||vS||^2>0. Actual deletion i gives Q<1.
The old affine witnesses exclude own3,4,6 and parallel contact just
as before, so gamma>=cos(pi/5), finite or dotted. Write

    c=A*p-gamma*omega_i-vS,
    d0=1-Q, u=d0/gamma^2,
    A^2=gamma^2*(D_i-u)>0.

The vector vS is a nonnegative combination of old foreign normals
by inverse-Gram positivity. If u<1, form

    xi=c+vS+(d0/gamma)*n_i.

Its selected products are xi.c=0, xi.foreign=0,
xi.n_i=gamma*(u-1)<0; its norm is -gamma^2*u*(1-u)<0.
All coefficients are positive, so every other actual wall inequality
also holds. Its positive A coefficient makes it future. Farthest
comparison and x.n_i>0 give the PURE B pressure

    h_i>sqrt(D_i-u)-sqrt(u*(1-u)).                    (1)

At u=1 the same vector is nonzero future null and satisfies every
wall inequality, giving an ideal point in the polytope's closure,
contrary to compactness.

If u>1, the actual finite set consisting of c, n_i and ALL foreign
old nodes has residual delta=gamma^2*(u-1)>0. The positive acute
projection gives
x_B.gamma*(omega_i-n_i)>=H(A,delta)>A-sqrt(delta). Hence

    h_i>sqrt(D_i-u)-sqrt(u-1).                       (2)

No projected mirror is assumed actual or native in either construction.
For D_i>=10, (1) has loss from sqrt(D_i) less than

    u/6+sqrt(u*(1-u)) <=(1+sqrt(37))/12<3/5.

For (2), Q>=1/4 gives u<=3/phi^2=6-3phi<7/6, with phi=(1+sqrt5)/2.
Thus sqrt(u-1)<=2-phi<2/5. Also
sqrt(D_i)-sqrt(D_i-u)=u/[sqrt(D_i)+sqrt(D_i-u)]<1/5:
the denominator is>3+2.9>35/6 and u<7/6. The total loss is<3/5.
Every foreign isolated cap is therefore strong, uniformly through
i=r-2. At Q0, the accepted isolated-cap arguments remain applicable.

## 3. Foreign far-end packets are strong

For far3@1 plus own gamma>=cos(pi/5), the full actual Schur data are

    delta=(i+1)/(2i)-Q>0,
    cosh L_i=(2i*gamma+1)/sqrt(i+1-2iQ).

For Q>0, Q>=1/4. The independently reviewed finite-index lemma
`compact_classical_foreign_far_packet_strong_20261020.md`, SHA256
`3766485d02300a78951e26e7a162086aa4ff74db9dc3aaa0842584ffb5fb1520`,
proves h_i>sqrt(2i)-.61 for ALL5<=i<=29. Its exact comparisons retain
Q and use the actual half-tangent identity; no asymptotic estimate is
substituted. Our indices satisfy this range because r<=29.

Consequently EVERY weak5<=i<=r-4 has exactly the old far-helper packet

    v=gamma*omega_i+omega1/2,
    A^2=2i*gamma^2+2gamma-1/2,
    Q=0, own label5,>=7 or dotted,
    101/125<gamma<6/5.                              (3)

## 4. The two short-tail layers, including zero-own branches

At i=r-3 the additional right B3 extension is short3@r making F4.
It seals every foreign/left contact. Own0 gives affine F4; own>=3
has delta=1/4,sigma>=1/2+1/(2sqrt2), hence is strong. All other
branches are the preceding ones, so (3) extends through r-3.

At i=r-2 the extra right B2 extension is short3@r. With no left
contact, own0 gives affine F4, own4 the affine4,3,3 star, own6 the6,3
path, and own>=5 is strong. The remaining own3 profile can have Q.
If Q>0, its delta=1/2-Q<=1/4 and
sigma=1/2+1/(2sqrt2), again strong. Thus its ONLY weak version is

    Q=0, v=(omega_(r-2)+omega_r)/2.                 (4)

If the extra short3 joins the left path as well, finite classification
forces a simple3 at a left endpoint and the full connected white B_r;
this seals every foreign contact. Its delta=1/(2i) and
sigma>=gamma+1/(2i)+1/(2sqrt2), strong for i>=21 even at gamma0.
Thus the zero-own short-B2 bridge is retained and estimated, not
discarded by the old whole-vertex positive-index argument.

At Q0, isolated caps at these layers still have the LEFT actual H4
projection with vector gamma*(3omega_(i+1)/4+omega_(i-3)/4).
Since i+1<=r-1 is long, prefix concavity bounds its x-product by
gamma*h_i. All other isolated/core-foot inputs also remain available.
Therefore (3)--(4) are exactly the accepted weak-tail inventory.

## 5. Actual-point and positive-dual closure

The complete weak inventory now has Q0. Every mixed actual vertex,
negative-core point and finite pressure used in the two accepted
whole-B sources therefore has displacement supported in B. Retaining
all old spectator walls adds a separate positive definite block.
Positive acute projection still makes each selected point ACTUAL
against every unknown wall. Prefix concavity and Cauchy need only
||x_B||<1, which follows from the global farthest vector.

Hence the accepted proofs apply with no changed constant: packets(3)
are nonadjacent; index5 is strong, using strong6 and one of16/17;
the profile(4) has its actual strong r-1 mate for r23..25; and a weak
far21 forces x1>.424. The latter actual H4/core feet need only simple
neighbors of i, available through r-2. At r23 the repaired energy is

    ||x_B||^2>.99752+(5/4)*(.424-.361)^2-7/3200
                =1.00029375>1.

All other r23..25 cases have the accepted endpoint21/22/23 reference
and its controlled isolated/two-hole loss. For r>=26 use the accepted
endpoint21/22 reference. All index choices are inside B; none uses a
spectator rate. This contradicts ||x_B||<1 and completes the theorem.

The generic actual-shadow/positive-point inputs are SHA256
`c913ec3bd8e0f26a8c679a1ede2b98b251ad3f36efeae60486213ddfa3e8badc`,
`3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`,
and orthogonal-component opposite obstruction SHA256
`4216d06c2e03fd336898b8dac289b9969fa096c8a1d08bd67ddeb8ab5ec2f23d`.
The general isolated construction is also written independently in
Section2 of `compact_an_large_components_arbitrary_spectator_transfer_20261020.md`,
SHA256 `8d3b999761be7dd46e7990973b3ec4d4a746d4d682b044d53695f4aae405e3e5`.
No arithmetic, literal-host exclusion or new ambient ceiling is used.
