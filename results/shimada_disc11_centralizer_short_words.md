# Short words in the enhanced discriminant-11 K3 centralizer

## Exact bounded result

Let \(X\) be Shimada's discriminant-11 singular K3 surface, let
\(\Gamma=\operatorname{im}(\operatorname{Aut}(X)\to O(\operatorname{NS}(X)))\),
and use Shimada's published generating alphabet

\[
\mathcal G=\{g_0,\ldots,g_{766}\}.
\]

For the Weyl-framed golden correspondence \(\tau\) whose eigenspace meets the
physical chamber \(D_{68}\), exhaustive exact search gives

\[
\begin{array}{c|c|c}
\text{word length}&\text{centralizing ordered words}&
\text{distinct actions}\\ \hline
1&0&0\\
2&682&\{1\}\\
3&4989&\{1,A\},
\end{array}
\]

where \(A\) is exactly the previously discovered D68/D27 return involution.
In Shimada's zero-based generator indexing it has the direct surface-group
factorization

\[
\boxed{A=g_{507}g_1g_{39}.}
\]

It satisfies

\[
A^2=1,qquad
\chi_A(t)=(t-1)^{18}(t+1)^2.
\]

Thus the bounded surface-group search finds no infinite-order enhanced
automorphism and no second finite generator through word length three.  This
is evidence about a completely specified generating alphabet, not a proof
that the full centralizer is \(C_2\).

## Exact meet-in-the-middle certificate

Shimada's matrices act on Neron--Severi row vectors.  Put
\(C_i=g_i^{\mathsf T}\).  A two-letter word centralizes \(\tau\) precisely
when

\[
C_jC_i\tau=\tau C_jC_i,
\]

or equivalently

\[
C_i\tau C_i^{-1}=C_j^{-1}\tau C_j.
\]

Hashing these two exact conjugacy lists exhausts all \(767^2\) ordered words.
The 682 matches all multiply to the identity.  Their ordered-pair digest is

```text
3c4bd2e1d1893c3f7c04aa00bec9850be30f367204661b9e56788a09be214eb5
```

For length three, the necessary equation is

\[
C_j(C_i\tau C_i^{-1})C_j^{-1}=C_k^{-1}\tau C_k.
\]

A modular two-vector probe over \(\mathbf F_{1000003}\) filters all
\(767^3\) triples by evaluating both sides.  This is only a necessary filter:
every surviving triple is reconstructed over \(\mathbb Z\), multiplied as a
full \(20\times20\) matrix, and its commutator with \(\tau\) is checked
entrywise.  Modular collisions can therefore add work but cannot create a
false answer or discard an equality.  The 4,989 accepted triples have digest

```text
af277e3506b6d0da7d9dfd5f076361b848be71e1a6a790bbe4cd40c2a08f0f47
```

and collapse to exactly the identity and \(A\).  The independently derived
D68 gallery factorization \((507,1,39)\) is required to occur in this exhaustive
set, providing a cross-check between the surface chamber atlas and the global
generator list.

## Relation to the algebro-geometric decision theorem

The enhanced covolume theorem proves

\[
\operatorname{vol}P_{X,\tau}<\infty
\iff P_{X,\tau}\text{ compact}
\iff |\operatorname{Aut}(X,\tau)|<\infty.
\]

Consequently a proof that the complete centralizer is \(\langle A\rangle\)
would certify a compact Coxeter \(9\)-polytope.  Conversely, one
infinite-order word would certify infinite volume.  The present word-length
bound decides neither alternative; its role is to calibrate the direct
algebro-geometric centralizer search and eliminate all shorter candidates.

## Reproduction

With the SHA-pinned `disc11.txt` supplement available, run

```sh
PYTHONPATH=src:. python3 \
  scripts/search_shimada_disc11_centralizer_words.py /path/to/disc11.txt
```

The compact certificate is
`data/candidates/shimada_disc11_centralizer_words3.json`.
