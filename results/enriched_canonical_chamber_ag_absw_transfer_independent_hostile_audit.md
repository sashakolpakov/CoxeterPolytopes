# Independent hostile audit: enriched canonical-chamber AG/ABSW transfer

**Verdict (2026-08-11): PASS CONDITIONAL after five scope repairs.**  The
dimension theorem is conditional on exhaustive entry into one of three
fully certified carrier branches.  The finite-field theorem additionally
assumes an exact full-Picard complex K3 trace carrier and a class-uniform
mixed correspondence height.  Neither raw essential surjectivity nor that
height implication has been proved in the non-rational \(n\ge2\) lane.

This audit does not import the owner certifier.  It independently pins the
branch theorems and their hostile audits, rebuilds every field constant, and
checks the owner archive only after the independent calculation.

## 1. The three dimension branches do not leak hypotheses

The direct rational row uses a literal full Neron--Severi nef chamber, the
complete negative-extremal \((-1),(-2),(-4)\) alphabet, integral facet
reflections, and host reflection maximality.  These are the hypotheses of
the pinned free-chain/MMP theorem, which gives \(n\le11\) and, because all
literal curve normals are rational, \(k=\mathbf Q\).

The weak-log row is now explicitly a **smooth projective complex rational**
carrier.  It requires a disjoint \((-4)\) boundary, actual irreducible
proposed facets, the adapted host, relative maximality, and an exact
basepoint-free certificate for

\[
                         D=-(K+B/2).
\]

The basepoint-free certificate makes \(D\) nef independently of the proposed
wall list.  The pinned exact-return and orbifold second-Chern theorem then
gives \(n\le25\) and \(k=\mathbf Q\).

The non-rational row returns the entire restriction-of-scalars block and
canonical-preserving trace operators on a minimal carrier, or supplies a
geometric contraction with finite numerical exceptional orbit.  The pinned
descent theorem places the full block in the canonical orthogonal of a
minimal Kodaira-dimension-zero complex surface.  Hence

\[
                  d(n+1)=\operatorname{rank}T\le20.
\]

This gives \(n\le9\) when \(d\ge2\).  Taking the largest rational row gives
the conditional global ceiling \(n\le25\).  The standard marked ordinary
lane has the separately audited exact property-4 theorem with \(n\le9\)
and explicit carriers for every \(3\le n\le9\).

## 2. The field theorem really is full-Picard K3

Gate 5 does not merely embed the trace block in a larger Picard lattice.  It
requires an exact full-Picard complex K3 carrier, so

\[
                 \operatorname{NS}(X)=T,
                 \qquad \rho(X)=d(n+1)\le20,
\]

and returns \(H\in T\) big and nef.  The trace form is

\[
                 q=\operatorname{Tr}_{k/\mathbf Q}(\alpha h),
\]

with one Lorentzian place and negative-definite non-distinguished
summands.  The cycle Cauchy--Schwarz argument uses the positive-definite
form \(h^\sigma\), equivalently \(-q_\sigma\); the repaired synthesis now
states this sign explicitly.

For a dominant totally positive \(\theta\), write the non-distinguished
contributions of \(H\) as \(-a_\sigma\), with \(a_\sigma\ge0\).  Then

\[
\begin{aligned}
 q(H,\tau_\theta H)-\sigma_0(\theta)q(H,H)
   &=\sum_{\sigma\ne\sigma_0}
     \bigl(\sigma_0(\theta)-\sigma(\theta)\bigr)a_\sigma\\
   &\ge0.
\end{aligned}
\]

Thus the assumed mixed bound and \(H^2\ge s_0\) give
\(\sigma_0(\theta)\le C/s_0\).  Definite-place Cauchy--Schwarz bounds the
other conjugates.  This is the exact point at which the scalar
correspondence height enters; Picard rank or \(H^2\) alone does not replace
it.

## 3. Cycle length, generators, and the discriminant exponent

Trace induction on a finite projective returned order-module makes all
directed cyclic Cartan products integral.  Pair cycles and chordless simple
cycles generate the field.  At a non-distinguished place, the chordless
cycle Gram kernel has dimension at most two, so

\[
                          \ell\le n+3.
\]

The shifted coordinate

\[
                \theta_C=(-1)^\ell c(C)+2^\ell
\]

is totally positive.  It is either bounded at all places by
\(2^{\ell+1}\), or dominant and controlled at the distinguished place by
\(C/s_0\).  Hence every selected coordinate has all conjugates in
\((0,R_{0,n}]\), where

\[
 R_{0,n}=\max\{2^{n+4},C/s_0\}.
\]

If \(D=\lfloor20/(n+1)\rfloor\ge2\), at most
\(r=\lfloor\log_2D\rfloor\) successive field generators are needed:
every strict step multiplies degree by at least two.  There are at most
\(D(D-1)/2\) embedding-pair hyperplanes on which a linear combination is
not primitive.  A positive grid of side

\[
                       M_D=D(D-1)/2+1
\]

cannot be covered by them, because each hyperplane contains at most
\(M_D^{r-1}\) grid points.  A primitive integral combination therefore has
all conjugates in

\[
                   (0,R],\qquad R=rM_DR_{0,n}.
\]

Positivity matters: each conjugate difference is at most \(R\), not
\(2R\).  The squared Vandermonde contains \(D(D-1)\) ordered difference
factors, so

\[
                  |\operatorname{Disc}k|
                  \le R^{D(D-1)}.
\]

The elementary symmetric functions then put every coefficient of the
monic minimal polynomial in a finite integer box.  This proves a literal
finite field set, not merely bounded degree.

## 4. The constants (8192), (32), and (30)

For \(2\le n\le9\), the independent table is

\[
\begin{array}{c|c|c|c|c|c}
n&D&n+3&rM_D&2^{n+4}&D(D-1)\\ \hline
2&6&5&32&64&30\\
3&5&6&22&128&20\\
4&4&7&14&256&12\\
5&3&8&4&512&6\\
6&2&9&2&1024&2\\
7&2&10&2&2048&2\\
8&2&11&2&4096&2\\
9&2&12&2&8192&2.
\end{array}
\]

Therefore every dimension-specific base is at most
\(32\max\{8192,C/s_0\}\), and every exponent is at most \(30\).  This
proves the coarse uniform formula

\[
 \left(32\max\{8192,C/s_0\}\right)^{30}.
\]

At \(n=2\), the exact row gives

\[
 \left(32\max\{64,C/s_0\}\right)^{30}.
\]

The repaired theorem first branches on \(D_n\).  When \(D_n<2\), a
non-rational full-trace K3 carrier is impossible, while the rational
branches already have \(k=\mathbf Q\).  In particular, for
\(20\le n\le25\), the K3 trace ceiling is zero but the actual field degree
is correctly recorded as one.  No \(\log_2(0)\) is formed.  The archive now
replays all dimensions \(2\) through \(25\).

## 5. Signed Hilbert/Chow grammar

A signed Chow class is not itself a point of one ordinary Chow variety.
The repaired Corollary 3.2 therefore requires each signed marked basis
class and each signed correspondence to be returned as a difference of
positive and negative effective parts.  It bounds the total degree of
**both** parts and every rational denominator.  Together with a finite
embedding-type list, a finite return grammar, and the requirement that the
marked basis span the exact Picard lattice on every accepted stratum, this
gives finitely many smooth-proper strata.  Intersection and action matrices
are locally constant on each stratum, so both the mixed height and the
finite-host clause follow.  A bounded polarization square or a per-component
degree bound would not suffice.

## 6. Fields, hosts, and chambers remain separate conclusions

The finite-field theorem ends before integral hosts are counted.  A field
discriminant bound does not bound the discriminant or isometry type of a
reflection lattice.  Chamber finiteness additionally assumes finitely many
host isometry classes or complete normalized packets, preservation of the
host, and reflection maximality relative to it.  Only then does
\(W(P)=R(M)\), and local finiteness plus connectedness of the complete
mirror chamber graph gives one chamber orbit per host.

The rank-four K3 counterfamily is in hyperbolic dimension one.  It proves
that the first four weak axioms do not control spectral height once
\(n=1\) is allowed, but it does **not** prove Gate 5 logically necessary in
the present \(n\ge2\) scope.  Whether exact maximal eigen-nef return forces
the required height in higher dimension remains open.  The other open
essential-surjectivity boundary is the rational twice-primitive \((-4)\)
lane.

The independent executable audit is
[audit_enriched_canonical_chamber_ag_absw_transfer_independent_hostile.py](../scripts/audit_enriched_canonical_chamber_ag_absw_transfer_independent_hostile.py),
with separate tests in
[test_enriched_canonical_chamber_ag_absw_transfer_independent_hostile.py](../tests/test_enriched_canonical_chamber_ag_absw_transfer_independent_hostile.py).
