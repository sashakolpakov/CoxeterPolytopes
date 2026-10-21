# Projected K3 scale and glue do not obstruct maximal abstract-chamber finiteness

**Status (2026-08-11): proved, carrier-relative, using Nikulin's
fixed-field/fixed-rank reflective-lattice finiteness theorem.**  Finite
ambient Neron--Severi hosts are not needed to count normalized abstract
chambers.  The correct intermediate object is the native arithmetic
reflection lattice **up to similarity**.  Global rescaling of its form does
not change a reflection, a normalized Coxeter Gram matrix, or a hyperbolic
chamber.

Let \(\mathfrak C_{\mathrm{proj}}\) be any class of finite-sided,
finite-volume arithmetic Coxeter chambers \(P\subset\mathbf H^n\), \(n\ge2\),
with the following properties.

1. Each \(P\) has an exact projected complex-K3 eigenslice return containing
   the full native trace module.  Thus, for its totally real ground field
   \(k\), degree \(d=[k:\mathbf Q]\), and native rank \(m=n+1\),
   \[
                         dm\le \rho(X)\le20.             \tag{0.1}
   \]
   A rational return may be a proper projector slice in a larger Picard
   lattice.  A nonrational return retains the full trace and admissibility
   data.
2. The facet group \(W(P)\) is maximal among discrete hyperbolic reflection
   groups.  Maximality among arithmetic reflection groups also suffices.

Then only finitely many unmarked normalized abstract chambers occur in
\(\mathfrak C_{\mathrm{proj}}\).  The ambient Picard lattices, orthogonal
complements, projector denominators, and glue indices may all be infinite in
number.

This gives the exact conditional ABSW composition:

\[
\boxed{
 \begin{gathered}
 \text{class-wide exact projected K3 return}\\
 +\ \text{raw reflection maximality}
 \end{gathered}
 \Longrightarrow
 \begin{gathered}
 n\le19,\quad\text{finite ground-field set},\\
 \text{finitely many normalized maximal chambers}.
 \end{gathered}}                                      \tag{0.2}
\]

The last implication is not a classification-free AG proof: it uses the
fixed-data theorem of Nikulin stated precisely in Section 3.  What is new in
the present composition is the sharp separation between ambient K3 hosts and
native similarity hosts.  The projected-root scaling family invalidates
finiteness of the former but not of the latter.

## 1. What projected K3 geometry supplies

Equation (0.1) immediately gives

\[
                         3\le m\le20,
             \qquad 2\le n=m-1\le19.                   \tag{1.1}
\]

For \(d>1\), it gives the sharper \(m\le10\), hence \(n\le9\).
The independently proved finite-volume K3 field theorem gives a finite
ground-field set:

- a noncompact exact admissible carrier has \(k=\mathbf Q\);
- a compact carrier has its field in the finite effective set cut out by
  \[
   |\operatorname {Disc}(k)|
       \le (32\cdot14^{10})^{30}.                       \tag{1.2}
  \]

No physical same-vector normalization is used here.  In particular, neither
(1.1) nor (1.2) bounds the scale of the trace form or the discriminant of the
ambient Picard lattice.

There is nevertheless a useful projected-root normalization.  Choose a
minimal-aperture closed star and extract \(m\) facets whose projected normals
span the native Lorentz space.  Let a physical facet root decompose
orthogonally as

\[
                 \delta_i=t_i+u_i,\qquad
 t_i\in T_{\mathbf R},\quad u_i\in T_{\mathbf R}^{\perp},
 \qquad \delta_i^2=-2.                                 \tag{1.3}
\]

Put \(a_i=-t_i^2\).  Visibility of the mirror gives \(a_i>0\), while the
negative definiteness of the complement gives

\[
                0<a_i\le2,\qquad -u_i^2=2-a_i.         \tag{1.4}
\]

Every two projected normals in the closed star satisfy

\[
            |t_i.t_j|<7\sqrt{a_i a_j}.                  \tag{1.5}
\]

Cauchy--Schwarz in the negative complement and (1.5) imply

\[
\begin{aligned}
 \delta_i.\delta_j
 &\le 7\sqrt{a_i a_j}
       +\sqrt{(2-a_i)(2-a_j)}\\
 &\le14.                                                \tag{1.6}
\end{aligned}
\]

For the last inequality, apply the operator
\(\operatorname {diag}(7,1)\) to the two vectors
\((\sqrt{a_i},\sqrt{2-a_i})\) and
\((\sqrt{a_j},\sqrt{2-a_j})\), both of Euclidean norm \(\sqrt2\).
The first inequality in (1.5) is strict, so (1.6) is strict.  Distinct
effective irreducible K3 facet curves have nonnegative integral
intersection.  Consequently the selected physical lifts satisfy

\[
                       0\le\delta_i.\delta_j\le13.       \tag{1.7}
\]

This lemma is genuinely scale-free, but it does **not** prove ambient-host
finiteness: the \(\delta_i\) need not span \(\operatorname {NS}_{\mathbf Q}\),
and their ordinary K3 root reflections need not be the raw native
\(k\)-reflections.  It is a bounded local certificate, not the global
chamber-counting step.

## 2. The raw arithmetic lattice is reflective

Fix one raw chamber \(P\), its admissible native \(k\)-space \((V,h)\), and
an integral \(\mathcal O_k\)-lattice \(L\subset V\) preserved by all facet
reflections.  Such a lattice is part of arithmeticity.  Starting with a
module over a suborder causes no problem: take its \(\mathcal O_k\)-span and
multiply the form by one common element of \(k^\times\) to clear the finite
denominators.  This changes neither the reflections nor the hyperbolic
space.

Let \(O^+(L)\) be the integral isometry group preserving the chosen positive
cone, and let

\[
 R(L)=\langle s_r:s_r(L)=L\text{ and }s_r
          \text{ is a native hyperplane reflection}\rangle.              \tag{2.1}
\]

Admissibility makes the nondistinguished real orthogonal factors compact, so
the distinguished image of \(O^+(L)\) is an arithmetic lattice.  We have

\[
                         W(P)\subseteq R(L)\subseteq O^+(L).              \tag{2.2}
\]

Both \(W(P)\) and \(O^+(L)\) have finite covolume.  Therefore
\([O^+(L):W(P)]<\infty\): a fundamental domain for a subgroup consists,
up to null boundaries, of one translate for every coset, so an infinite
index would give infinite covolume.  It follows from (2.2) that

\[
                         [O^+(L):R(L)]<\infty.           \tag{2.3}
\]

Thus \(L\) is a reflective hyperbolic lattice in the exact sense needed by
the fixed-data theorem.  Notice that no K3 root need lie in \(L\).

## 3. Fixed field and rank leave finitely many native similarity hosts

The arithmetic input is the following theorem.

### Theorem 3.1 (Nikulin fixed-data reflective-lattice finiteness)

For a fixed totally real field \(k\) and fixed rank \(m\ge3\), there are only
finitely many similarity classes of hyperbolic \(\mathcal O_k\)-lattices
\(L\) for which \([O(L):R(L)]<\infty\).

This is Theorem 5.2.1, pp. 119--122, of V. V. Nikulin, *On the
classification of arithmetic groups generated by reflections in Lobachevsky
spaces*, Math. USSR Izv. **18** (1982), 99--123.  The pinned primary PDF has
SHA-256
`696652b6a67c6ed420de6c1cb61d9422e3619e28e081a78131a5323ce2edee9d`.
The project's [conceptual source audit](nikulin_conceptual_dependency_hostile_audit.md)
separately checks the theorem's role from the later dimension and field
arguments.

For clarity, here is the proof mechanism rather than a black-box slogan.
The geometric appendix used by Nikulin selects roots
\(r_1,\ldots,r_m\in L\) such that

1. \(\Delta=\sum_i\mathcal O_k r_i\) has finite index in \(L\);
2. the nonzero-pairing graph of the \(r_i\) is connected; and
3. for one absolute constant \(A\) (the source takes \(A=62\)),
   \[
      c_{ij}={4(r_i.r_j)^2\over r_i^2r_j^2}<A^2          \tag{3.1}
   \]
   at the distinguished place.

Reflection integrality makes every \(c_{ij}\) an algebraic integer.  At a
nondistinguished definite embedding, Cauchy--Schwarz gives
\(0\le\sigma(c_{ij})\le4\).  For fixed \(k\), the Minkowski image of
\(\mathcal O_k\) is a lattice, so only finitely many \(c_{ij}\) can occur in
these bounded boxes.  Hence there are finitely many connected normalized
types \((c_{ij})\).

Write \(d_{ij}=r_i.r_j\).  The root conditions are

\[
                         d_{ii}\mid 2d_{ij}.             \tag{3.2}
\]

For a prime ideal not dividing \(2\) or a nonzero \(c_{ij}\), (3.1)--(3.2)
force the valuations of \(d_{ii},d_{ij},d_{jj}\) to agree along every
nonzero edge.  Connectedness therefore makes the valuation common to the
whole matrix; it is removed by one global similarity.  At the finitely many
bad primes, the same equations bound all valuation differences by a
constant depending only on \(m\), \(v(2)\), and the finite type
\((c_{ij})\).  Thus only finitely many divisor matrices occur modulo a
common divisor.

After the divisor matrix is fixed, the remaining unit matrices are finite
modulo rescaling the individual \(r_i\): the square-unit subgroup has finite
index in \(\mathcal O_k^\times\), and changing
\(r_i\mapsto\varepsilon_i r_i\) absorbs those squares.  This gives finitely
many similarity classes of the full-rank frame lattice \(\Delta\).
Finally, reflection integrality gives the uniform index estimate printed in
the source,

\[
 [L:\Delta]\le [\Delta^*:\Delta]\,
                   [\mathcal O_k:2\mathcal O_k]^m,       \tag{3.3}
\]

after choosing a representative of the similarity class.  Hence only
finitely many intermediate projective lattices \(L\) occur.  This proves
Theorem 3.1.

The connected frame in this proof is the missing global object.  The
minimal-aperture projected K3 star gives a bounded spanning packet, but its
raw Gram graph can have orthogonal components (a right-angled facet already
shows why).  Concretely, in a regular compact right-angled hyperbolic
pentagon every side has minimal aperture.  The closed star of one side is
that side together with its two neighbors.  Its three normals span the
Lorentz space, but the chosen side is orthogonal to both neighbors, so the
nonzero-pairing graph is disconnected.  Thus one cannot silently replace
Theorem 3.1 by the local star bound; the relative scale between its
orthogonal blocks is not controlled by that packet.

## 4. Maximality collapses every similarity host to one chamber

Apply Theorem 3.1 to the reflective lattice obtained in Section 2.  For each
fixed pair \((k,m)\), only finitely many native similarity hosts \([L]\)
occur.  Similarity is exactly the right equivalence: replacing \(h\) by
\(a h\) gives

\[
 x-{2a h(x,r)\over a h(r,r)}r
      =x-{2h(x,r)\over h(r,r)}r                         \tag{4.1}
\]

and leaves every normalized Gram entry unchanged.

Now use raw reflection maximality.  The complete group \(R(L)\) is a
discrete reflection overgroup of \(W(P)\).  Hence

\[
                              W(P)=R(L).                 \tag{4.2}

\]

If maximality was assumed only among arithmetic reflection groups, (4.2)
still follows: by (2.3), \(R(L)\) is an arithmetic lattice commensurable with
\(W(P)\).

The complete mirror arrangement of \(R(L)\) is locally finite.  A generic
compact path between regular points in two chambers crosses finitely many
mirrors; reflecting successively in those mirrors carries one chamber to
the other.  Thus the complete arrangement has one \(R(L)\)-orbit of
chambers.  Each oriented native similarity host contributes at most one
unmarked normalized chamber class.

There are finitely many fields by Section 1 and finitely many ranks by
(1.1).  Taking the finite union over \((k,m)\), Theorem 3.1 and (4.2) proves
the assertion in (0.2).

### Corollary 4.1 (finite-pair principle)

Let \(\mathfrak C\) be **any** class of finite-volume arithmetic Coxeter
chambers of dimension at least two, maximal among discrete or arithmetic
reflection groups.  If

\[
              \{(k(P),\dim P):P\in\mathfrak C\}          \tag{4.3}
\]

is finite, then \(\mathfrak C\) contains only finitely many unmarked
normalized chamber classes.

Indeed, Sections 2--4 use the carrier only to obtain the finite set (4.3).
Once a pair \((k,m)\) is fixed, Theorem 3.1 finite-izes the native similarity
hosts and maximality gives (4.2).  Consequently an exhaustive revised AG
carrier stack would recover full ABSW finiteness as soon as its branches
give a finite dimension set and a finite ground-field set.  The existing
direct and weak-log rational branches have \(k=\mathbf Q\); exact
finite-volume full-trace K3 branches in \(n\ge2\) have the proved finite
field set.  What remains upstream is class-wide exact return into a revised
exhaustive stack, not finite ambient hosts.

## 5. The two hostile controls

### 5.1 Ambient host finiteness is false but irrelevant

The compact golden \((2,5,5)\) family has trace form \(sQ\), glue index
\(s^3\), and ambient discriminant \(64s^3(s-1)^3\) for every \(s\ge2\).
Its physical walls lie outside the trace block.  Nevertheless,

\[
                        [sQ]=[Q]                         \tag{5.1}
\]

as native similarity hosts, and every normalized Cartan product and the
abstract triangle are fixed.  This is exactly why the
[projected-root scaling no-go](k3_compact_projected_root_scaling_no_go.md)
does not contradict the present theorem.

### 5.2 Reflection maximality is essential

Fixed rank, field, and even a fixed divisor lattice do not bound arbitrary
finite-index reflection-subgroup chambers.  The
[K3 \(H^{10}\) redoubling family](k3_s3_h10_redoubling_family.md) has, for
every \(N\ge1\), a finite-volume arithmetic Coxeter chamber \(P_N\) on one
fixed rational rank-eleven lattice, with

\[
 \#\operatorname {Facets}(P_N)=N+12,
 \qquad [W(P_1):W(P_N)]=N,
 \qquad \operatorname {vol}(P_N)=N\operatorname {vol}(P_1).             \tag{5.2}
\]

The facet counts make the normalized abstract chambers pairwise
nonisometric as \(N\) varies.  These are descendants, not reflection-maximal
groups.  Equation (5.2) proves that the maximality step in (4.2) is doing
real work; no theorem based only on bounded \((k,n)\) can omit it.

## 6. Exact ABSW and AG-only scope

The proved conclusion is:

> **Carrier-relative maximal chamber finiteness.**  Exact projected K3
> full-trace returns in dimension at least two, together with raw reflection
> maximality, have only finitely many normalized abstract chambers.  Finite
> ambient Picard hosts are unnecessary.

If a broader projected-K3 compiler were exhaustive for all raw maximal
arithmetic chambers in its declared domain, this theorem would turn its K3
rank and field conclusions into ABSW-type finiteness on that domain.  A
K3-only compiler cannot cover dimensions above \(19\), and its class-wide
existence is independently open.

This is not yet an AG-only derivation of ABSW.  Remove Theorem 3.1 and the
present project has only:

\[
 \text{finite }(k,n)
 +\text{bounded local projected/Cartan star}
 \quad\not\Rightarrow_{\text{proved here}}\quad
 \text{finite complete raw similarity hosts}.           \tag{6.1}
\]

An independent AG replacement must produce the connected full-rank raw
reflection frame of (3.1), a finite native projective order, a bounded
orbifold volume/facet invariant, or an equivalent global completion datum.
The bounded physical intersections (1.7) alone do not suffice because the
physical root reflections are not the raw reflections and the roots need
not span the ambient Picard lattice.

## 7. Reproduction

- [certificate generator](../scripts/certify_projected_k3_maximal_abstract_chamber_finiteness.py)
- [exact certificate](../data/enumerations/projected_k3_maximal_abstract_chamber_finiteness.json)
- [focused regressions](../tests/test_projected_k3_maximal_abstract_chamber_finiteness.py)

Run

```text
python3 -m scripts.certify_projected_k3_maximal_abstract_chamber_finiteness --check-output
python3 -m pytest -q tests/test_projected_k3_maximal_abstract_chamber_finiteness.py
```
