# Degree-29 reduction for the `m=5` Gram--Belyi packet

**Status (2026-08-10): exact finite reduction; sharpness remains open.**

**Follow-up.**  The `1,372` three-black/`WW` records left open by this first
archive have since been eliminated exactly.  Thus any degree-`29` witness
must have at least four black vertices; see the
[three-black/`WW` no-go](gram_belyi_m5_degree29_three_black_ww.md).  The body
below retains the first-stage counts so that its frozen archive remains
auditable.
For

\[
 \lambda=\frac{7+2\sqrt5}{29},\qquad
 29\lambda^2-14\lambda+1=0,
\]

Rodriguez's theorem gives `h_Bel(lambda)>=29`.  The present calculation does
not find a degree-`29` normalized Belyi polynomial and does not prove that none
exists.  It does prove the following necessary condition.

> **Low-vertex reduction.**  Suppose a normalized degree-`29` Belyi
> polynomial maps `lambda` to zero or one.  Swap `B` and `1-B` so that
> `B(lambda)=0`, and call that fibre black.  The black fibre has at least
> three vertices.  If it has exactly three vertices, both normalized source
> points `x=0,1` are white.

It also excludes, for every integral exponent choice, the natural
Galois-symmetric family

\[
 B(x)=C(29x^2-14x+1)^a x^b(x-1)^c,\qquad
 2a+b+c=29.
\tag{1}
\]

This is a substantial exact ansatz no-go, not a no-go for all degree-`29`
dessins.

## 1. Passport census before coefficient solving

A degree-`29` Shabat polynomial has a `29`-cycle above infinity.  If its
black and white ramification partitions have `k` and `l` parts, the genus-zero
cycle equation is

\[
                         k+l=30.
\tag{2}
\]

Enumerating partitions of `29` subject to (2) gives `235,465` ordinary
Riemann--Hurwitz-compatible black/white passport pairs.  We then distinguish:

- the black cycle containing `lambda`;
- the colors and ramification indices of the ordered source points `0,1`;
- but not identities of cycles having the same ramification index.

This gives exactly `42,610,678` marked ramification passports.  These are
passport records, not `42,610,678` dessins and certainly not that many
polytopes.

For the first three black-vertex counts the exact ledger is:

| black vertices `k` | `0,1=BB` | `BW` | `WB` | `WW` | total | status |
|---:|---:|---:|---:|---:|---:|---|
| 1 | 0 | 0 | 0 | 1 | 1 | eliminated |
| 2 | 0 | 56 | 56 | 84 | 196 | all eliminated |
| 3 | 756 | 1,512 | 1,512 | 1,372 | 5,152 | `BB/BW/WB` eliminated here; `WW` eliminated in follow-up |

Thus `3,977` marked passport records are eliminated and `42,606,701`
remain.  Of the remainder, `1,372` are the three-black/`WW` stratum and
`42,605,329` have at least four black vertices.  The full histogram for
`k=1,...,29` is frozen in the machine archive.

The size of this remainder explains why an indiscriminate plane-tree search
would be the wrong next move.  The next bounded target is the `1,372`-record
three-black/`WW` coefficient system, whose opposite-color ramification
partition is necessarily either

```text
(3,1^26)       or       (2,2,1^25).
```

## 2. Exact low-vertex elimination

### One black vertex

There is only one black root, so

\[
                         B(x)=C(x-\lambda)^{29}.
\]

Both `0` and `1` must map to the other Belyi value.  Their quotient would
force

\[
                 \left(\frac{1-\lambda}{-\lambda}\right)^{29}=1.
\]

The ratio is negative real, while its odd power would have to be positive.
This eliminates the single marked passport.

### Two black vertices

Write

\[
                         F(x)=(x-\lambda)^a(x-t)^b,
                         \qquad a+b=29.
\tag{3}
\]

Its only non-root critical point is

\[
                         z=\frac{at+b\lambda}{29}.
\]

If one endpoint is the second black root, the other endpoint must satisfy
`F(endpoint)=F(z)`.  The certificate checks all `28` multiplicity choices on
both sides and finds a nonzero exact obstruction.

If both endpoints are white, one must have

\[
                         F(0)=F(1)=F(z).
\tag{4}
\]

For each `a=1,...,28`, (4) gives two polynomials in `t`, of exact degrees
`b=29-a` and `29`.  At a displayed good split prime their reductions have
the same degrees and gcd one.  Their characteristic-zero resultant is
therefore nonzero.  This eliminates all two-black passports, including the
possibilities that one white endpoint is the ramified white vertex.

### Three black vertices with a black endpoint

For a black endpoint at zero, write

\[
 F(x)=(x-\lambda)^a x^b(x-t)^c,qquad a+b+c=29.
\tag{5}
\]

The non-root factor of `F'` is

\[
 D(x)=29x^2-\bigl((a+b)t+(b+c)\lambda\bigr)x+b\lambda t.
\tag{6}
\]

All non-root critical points have a common value precisely when the remainder
of `F` modulo `D` is constant.  Write that remainder as `R_1(t)x+R_0(t)`.
The remaining white endpoint additionally requires

\[
 R_1(t)=0,qquad (1-\lambda)^a(1-t)^c-R_0(t)=0.
\tag{7}
\]

Both polynomials in (7) have exact degree `28`.  Here is the degree argument,
which matters for the modular certificate.  Homogeneous division of the
degree-`29` product in (5) by the homogeneous quadratic (6) makes `R_1`
homogeneous of degree `28`; the constant remainder is divisible by
`lambda`, since it vanishes after setting `lambda=0`, so its `t`-degree is
also at most `28`.  A good reduction has degree `28` in both coordinates,
proving that the characteristic-zero leading coefficients do not vanish.

For every one of the `378` triples `(a,b,c)`, a prime among
`101,109,131` gives gcd one in (7).  Applying `x -> 1-x` treats the other
endpoint and eliminates another `378` triples.  If both endpoints are black,
there is no parameter: the exact remainder coefficient is nonzero for all
`378` triples.  This proves the low-vertex reduction stated above.

## 3. Why the modular certificates are exact

The calculation is not using failure to find a point over a finite field.
For the one-parameter strata it does the following:

1. derive the two necessary polynomials over `Q(sqrt(5))`;
2. prove their characteristic-zero degree;
3. reduce at a split prime `p != 29` without losing either degree;
4. compute their full polynomial gcd over `F_p` and obtain `1`.

Gcd one is equivalent to a nonzero Sylvester resultant at that prime.  A
zero resultant over `Q(sqrt(5))` would reduce to zero at every good prime, so
one nonzero reduction proves that the original equations have no common
root over the algebraic closure.  This detects roots over `overline(F_p)`,
not merely `F_p`-rational roots.

The witness-prime histograms are:

| exact stratum | trials | prime histogram |
|---|---:|---|
| two black, one endpoint black | 56 | `101:54, 109:2` |
| two black, both endpoints white | 28 | `101:21, 109:7` |
| three black, both endpoints black | 378 | `101:367, 109:11` |
| three black, exactly one endpoint black | 756 | `101:736, 109:18, 131:2` |
| Galois-symmetric family (1) | 210 | `101:210` |

Every trial and its separating prime are retained in the archive; the table
is not a probabilistic sample.

## 4. Galois-symmetric quadratic-factor family

For (1), divide the derivative by all forced root factors.  The residual
logarithmic-derivative numerator has degree two, three, or less according to
whether the endpoints are themselves roots.  The polynomial is Shabat only
if all roots of this residual factor have one common nonzero value; equivalently
the remainder of the unscaled product modulo the residual factor is constant.
If `b=0` or `c=0`, the corresponding normalized endpoint must also have that
common value.

There are exactly

\[
                 \sum_{a=1}^{14}(30-2a)=210
\]

integer triples `(a,b,c)`.  Exact modular identity witnesses eliminate all
`210`.  This includes the conjugate quadratic root automatically and is
strictly stronger than testing only the distinguished real root inside this
ansatz.  It still says nothing about a general coefficient field in which
only `lambda`, not its conjugate, is forced into the black fibre.

## 5. Reproducibility and honest next step

The exact producer and archive are:

- [`certify_gram_belyi_m5_degree29_reduction.py`](../scripts/certify_gram_belyi_m5_degree29_reduction.py);
- [`gram_belyi_m5_degree29_reduction.json`](../data/enumerations/gram_belyi_m5_degree29_reduction.json).

One replay takes about one second and never launches a plane-tree census or a
numerical coefficient solve.  The archive explicitly records:

```text
degree_29_global_status = OPEN
degree_29_witness_found = false
complete_degree_29_no_go_claimed = false
```

The follow-up exact elimination has closed the three-black/`WW` stratum and
moved the search boundary to at least four black vertices.  Degree-`29`
sharpness itself remains open.
