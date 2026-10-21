# Independent hostile audit: projected-K3 rank saturation and the `X3` no-go

## Verdict

**PASS, with the owner's stated exact-module and non-maximal-only scope.**

The terminal owner package is internally correct on all five pressure
points:

1. a generic exposed K3 nef facet must be supported by an ambient
   square-\((-2)\) root;
2. rank saturation for the \(H^{19}\) `X3` chamber makes the projector the
   identity;
3. the pinned exact module permits only one rational common similarity, so
   the raw norm classes \(-2\) and \(-6\) would require \(3\) to be a rational
   square;
4. the compact nonrational residue shells are finite, and the rational
   cusped host-complete argument really does make every visible shell vector
   an old integral reflection; and
5. the facet-compression discriminant bound follows from primitive
   overlattice glue, not from facet count.

Accordingly `X3` disproves the statement

```text
every finite-volume arithmetic Coxeter chamber of trace rank <= 20
has an exact-module projected complex-K3 eigen-nef return.
```

It does **not** disprove a conjecture restricted to absolutely maximal
arithmetic reflection groups: the pinned `X3` chamber is complete for all
reflections of its discriminant-three host, but absolute maximality among
commensurable hosts is not proved.  Nor does it obstruct the existing
mixed-divisor K3 realization or an unlabelled-real carrier which forgets the
rational module.

Independent replay files are

- [audit script](../scripts/audit_projected_k3_universality_rank_saturation_no_go_independent_hostile.py),
- [audit archive](../data/enumerations/projected_k3_universality_rank_saturation_no_go_independent_hostile_audit.json),
  and
- [audit tests](../tests/test_projected_k3_universality_rank_saturation_no_go_independent_hostile_audit.py).

## 1. Frozen inputs and independent numerical replay

The owner files audited here have SHA-256 values

```text
7ed5176c0c0a36e036c0e12dfa6acf1ee15d90d5ec7143ca5df59ebc3b7143d2  report
d94e18eada49d277b811299f5378c6aa58566a861032c2d2d6fa78fdf3844c87  script
61ddd537135db3f47bd267abcc3d9bbea8db64d73370e7f8a1de2a6125e1751e  archive
9dcbd2189e5629bd8f251d9e0950c2e2a3d438f6a810f40e27bb56724b88ada1  tests
```

Without importing the producer, the audit reads the frozen `X3` wall Gram
and independently recomputes

\[
 \operatorname {rank}G=20,\qquad
 \operatorname {inertia}G=(1,19,16),\qquad
 \operatorname {SNF}_{\ne0}G=1^{19},3,
\]

together with the diagonal census

\[
                   24\times(-2)+12\times(-6).
\]

The audit also verifies the projector, self-adjointness, Schur complement,
root identity and residue-shell row bound on a separate rank-two,
three-facet block.  This is a formula audit, not a second proof of the frozen
Vinberg finite-volume enumeration.

## 2. Why a generic nef facet forces a root on its line

Let \(C\) be an ambient K3 root chamber and \(E\) the proposed eigenspace.
Take a positive-norm point in the relative interior of an exposed facet of
\(\mathbf P(C\cap E)\), away from all lower-dimensional faces.  K3 root
hyperplanes are locally finite in the positive cone.  In a small
neighborhood only finitely many root inequalities occur, so at least one of
them is active on an open subset of the facet.  Its projected normal is
therefore nonzero and lies on the facet line.

This rules out a possible escape in which infinitely many root mirrors
approach the facet but no individual root supports it.  The argument uses a
generic positive point, not an ideal endpoint; every hyperbolic facet has
such points even when the chamber is noncompact.

For `X3`, exact return has eigenspace dimension 20.  Since a complex K3 has
\(\rho\le20\), one gets

\[
                        \dim E=\rho=20.
\]

Thus \(E=\operatorname {NS}(X)_{\mathbb R}\), the orthogonal projector is
the identity, and each supporting ambient root itself lies on the raw facet
line.  There is no complement direction in which its norm could be
corrected.

## 3. The rational square-class contradiction

The exact-module contract retains the pinned rational lattice up to one
common rational similarity \(c>0\).  A rational K3 root on a raw rational
facet line is a rational multiple of the primitive raw normal.

Choose one raw square-\((-2)\) normal \(r\) and one raw square-\((-6)\)
normal \(u\).  If K3 roots on their lines are \(ar\) and \(bu\), with
\(a,b\in\mathbf Q^\times\), their physical squares give

\[
                    ca^2=1,\qquad 3cb^2=1.
\]

Hence

\[
                             (a/b)^2=3,
\]

impossible over \(\mathbf Q\), equivalently because \(v_3(3)=1\) is odd.
This uses exactly the common-\(\mathbf Q\)-similarity promise.  It would not
be valid for an unlabelled real polytope whose rational module had been
discarded, and the owner explicitly preserves that escape clause.

## 4. Compact nonrational residue shells

For one facet \(v_i\), the direct block reduces every visible root to a
vector \(y\) satisfying a fixed trace norm and one congruence:

\[
 B(y,y)=B(v_i,v_i),\qquad y\equiv\pm v_i\pmod{sT_0}.
\]

The owner's finiteness argument works over an arbitrary totally real field.
If the distinguished mirror meets compact \(P\), its distinguished normal
is spacelike.  Its contribution to the trace norm is therefore negative;
all nondistinguished contributions are negative because admissibility makes
those forms negative definite.  The fixed negative trace budget bounds the
nondistinguished components.  Compactness of \(P\) uniformly bounds the
Lorentz boosts of distinguished spacelike normals whose mirrors meet \(P\).
Thus the full Minkowski vector lies in a bounded set.  The Minkowski image of
the order module is a lattice, so only finitely many \(y\)'s occur.

Choosing \(s\) larger than every nonzero coordinate of the finitely many
\(y\mp v_i\) leaves only \(y=\pm v_i\) on the residue shells.  There is no
hidden assumption that equality of trace norms implies equality at every
field embedding.

## 5. The rational cusped host-completeness argument

For \(F=\mathbf Q\), let \(y=\pm v_i+sx\) be a residue-shell vector with
\(G(y,y)=G(v_i,v_i)\).  If \(s\) is divisible by
\(\lvert G(v_i,v_i)\rvert\), then for every \(z\in T_0\),

\[
 {2G(y,z)\over G(y,y)}
 =\pm {2G(v_i,z)\over G(v_i,v_i)}
   +{2sG(x,z)\over G(v_i,v_i)}\in\mathbf Z.
\]

The first term is integral because \(v_i\) is an integral reflection vector;
the second is integral by the divisibility choice.  Hence \(y\) is itself an
integral reflection vector of the pinned host.  If \(P\) is a chamber for
**all** such reflections, \(y^\perp\) cannot cut its interior.  This proves
residue-shell completeness even at a cusp.  The argument correctly does not
claim a nonrational analogue: equality of trace norms need not fix the
field-valued norm.

## 6. Facet compression and primitive glue

For a primitive trace lattice \(T\subset S\), its integral orthogonal
complement \(K=S\cap T_{\mathbb Q}^\perp\) is primitive as well.  Therefore
the overlattice glue

\[
                S/(T\oplus K)\subset A_T\oplus A_K
\]

projects injectively to both discriminant groups: it is the graph of an
anti-isometry between subgroups.  Each facet root decomposition
\(\delta_i=t_i+u_i\) contributes the glue pair
\(([t_i],[u_i])\).  The subgroup generated by the \([t_i]\)'s consequently
injects into \(A_K\).  Since the \(p\)-length of a rank-\(c\) lattice's
discriminant group is at most \(c=\operatorname {rk}K\),

\[
        c\ge\ell_p\bigl(\langle[t_1],\ldots,[t_m]\rangle\bigr).
\]

For uniform representatives \(t_i=v_i/s\) and \(p\mid s\), reduction of
their numerator classes maps onto the column span of \(V\bmod p\).  Hence

\[
                     c\ge\operatorname {rk}_{\mathbf F_p}(V\bmod p).
\]

This is a denominator/glue lower bound, not the false statement that the
number of facets lower-bounds complement rank.  The owner's compression
scope is correct.

## 7. Exact scope

The `X3` result is stronger than a failure of one proposed lift list and
weaker than a maximal-only Property-4 counterexample.

- It is stronger because local finiteness forces roots on the facet lines
  even when no preferred root witnesses are supplied.
- It disproves a universal theorem quantified over **all** finite-volume
  arithmetic chambers with exact trace rank at most 20.
- It does not address a theorem whose inputs are required to be absolutely
  reflection-maximal; that premise has not been proved for `X3`.
- It does not exclude the known `X3` mixed-divisor K3 carrier: its 24
  square-\((-2)\) facets are root walls, while its 12 square-\((-6)\) facets
  are extra integral divisor cuts.
- It gives no finite field or projected-host theorem.  The block theorems
  still retain `RSC(s)` and `EMB` outside their proved compact and rational
  host-complete subcases.

No scope repair to the owner bytes is required.

## Replay

```bash
python3 scripts/audit_projected_k3_universality_rank_saturation_no_go_independent_hostile.py \
  --check-archive

python3 -m pytest -q \
  tests/test_projected_k3_universality_rank_saturation_no_go_independent_hostile_audit.py
```
