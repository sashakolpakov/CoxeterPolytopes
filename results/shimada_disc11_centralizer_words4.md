# Length-four words in the enhanced discriminant-11 K3 centralizer

## Exact bounded theorem

Let \(X\) be Shimada's discriminant-\(11\) singular K3 surface, let

\[
\Gamma=\operatorname{im}\bigl(\operatorname{Aut}(X)\to
O(\operatorname{NS}(X))\bigr),
\]

and use Shimada's published alphabet
\(\mathcal G=\{g_0,\ldots,g_{766}\}\).  For the Weyl-framed golden
correspondence \(\tau\), exhaustive exact search proves:

\[
\boxed{\text{Every centralizing action represented by a word of length at
most four is either }1\text{ or }A.}
\]

Here \(A=g_{507}g_1g_{39}\) is the previously known D68/D27 return
involution.  In particular, no infinite-order enhanced automorphism occurs
through this bound.

The exact-length-four multiplicities are:

\[
\begin{array}{c|r}
\text{action}&\text{ordered syntactic words of length four}\\ \hline
1&1{,}138{,}938\\
A&439\\ \hline
\text{total}&1{,}139{,}377.
\end{array}
\]

These are word multiplicities, not distinct group elements.  A word of
syntactic length four may freely reduce or represent an action already
having a shorter expression.  The two matrices above are the only distinct
actions among all \(1{,}139{,}377\) accepted words.

## Exact meet-in-the-middle proof

Write \(C_i=g_i^{\mathsf T}\) for the column action.  For an ordered word
\((i,j,k,\ell)\), put

\[
A=C_jC_i,\qquad B=C_\ell C_k.
\]

Then the full column word is \(BA\), and

\[
BA\tau=\tau BA
\quad\Longleftrightarrow\quad
A\tau A^{-1}=B^{-1}\tau B.                 \tag{1}
\]

Thus the \(767^4=346{,}083{,}947{,}521\) ordered words are exhausted by
joining two lists of only

\[
767^2=588{,}289
\]

exact conjugates.  The implementation indexes the right-hand list in (1)
using the images of four fixed probe vectors modulo
\(p=1{,}000{,}003\), compressed by BLAKE2b-128.  This index is only a
necessary filter.  If (1) holds over the integers, the modular images and
their digest must agree, so a true solution cannot be omitted.  Digest or
modular collisions can only create extra candidates.

Every one of the \(1{,}139{,}377\) indexed candidates was reconstructed as
a full \(20\times20\) matrix over \(\mathbb Z\) using arbitrary-precision
FLINT arithmetic.  Its commutator with \(\tau\) was then checked entrywise.
All candidates passed; the modular stage was never used as evidence of
equality.  The accepted lexicographic word list has canonical compact-JSON
digest

```text
012f1d1e0314b1a70925cafa7e5e633971d86d55ff27d51393fe92b675df352c
```

Exact action hashing collapses the accepted words to matrices with digests

```text
identity: 8817db0b345cc277b1b2e649062c98d9a923e6a5159c08ff8b075a3947d0d259
A:        6b4bc88ee93385bcb4a4c55828924860aa90158f84a8469ffc238dbe4f33e4eb
```

The second is also compared entrywise—not merely by digest—with the
independently obtained length-three D68/D27 action, and satisfies

\[
A^2=1,\qquad \chi_A(t)=(t-1)^{18}(t+1)^2.
\]

## Scope

Together with the SHA-pinned length-at-most-three certificate, this is an
exact bounded theorem through word length four in the specified 767-letter
alphabet.  It does **not** prove that the full enhanced automorphism
centralizer is \(\langle A\rangle\): a new finite- or infinite-order action
could first occur at greater word length.  The global decorated-chamber
traversal remains the decision procedure for that unbounded question.

## Reproduction

With Shimada's SHA-pinned `disc11.txt` supplement and optional
`python-flint` available, run

```sh
PYTHONPATH=src:. python3 \
  scripts/search_shimada_disc11_centralizer_words4.py \
  /path/to/disc11.txt --progress
```

The machine-readable certificate is
`data/candidates/shimada_disc11_centralizer_words4.json`.
