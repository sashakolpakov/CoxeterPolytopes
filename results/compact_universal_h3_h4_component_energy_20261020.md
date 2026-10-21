# Short universal H3 and H4 component bounds

20 October 2026. Own written proof, awaiting outside-author review.
No scientific computation, optimizer, TeX build or Lean run was used.

At a vertex farthest from a strict interior point of a compact
hyperbolic Coxeter polytope, an H3 component has energy greater than
1/4 and an H4 component greater than 11/36. Other finite components
are arbitrary. These are component estimates, not dimension bounds.

Use the actual-wall and additive energy conventions of
`results/compact_universal_small_a_component_energy_20261020.md`,
SHA256 801581c4b855289aca7982b2d0f4bb6f5849f151a78f40b86c9d403e84f415e6.
In particular, with inverse old unit Gram C, an outside wall has height
square h, contacts beta, own cost q=beta^T C beta and directional
derivatives v_i=(C beta)_i/sqrt(C_ii). Then h=q+Q-1, where Q is the
sum of spectator costs. Each positive spectator cost is at least 1/4
and is strictly less than h: otherwise the same actual wall closes an
edge in that other component, contrary to the orthogonal-opposite
theorem. Equality would give a forbidden ideal edge endpoint.

Write phi=(1+sqrt(5))/2. All root lengths below are unit lengths.

## A common opposite under small energy

If the component energy is at most e, actual half-tangents satisfy
tau_i+tau_j <= sqrt(2e(1+k)), where k is their direction correlation.
The concavity of log(2t/(1+t^2)) gives

    r_i r_j <= [2sqrt(e(1+k)/2)/(1+e(1+k)/2)]^2,

where r_i=tanh L_i. Distinct actual opposites require r_i r_j>=k.
At e=1/4 the upper bound is 32(1+k)/(9+k)^2, strictly less than k
for every k>=1/sqrt(3). At e=1/3 it is 24(1+k)/(7+k)^2, strictly
less than k for every k>=1/sqrt(2): the cleared difference is
k^3+14k^2+25k-24, increasing, and its value at 1/sqrt(2) is positive
because 51/sqrt(2)>34.

The adjacent dual correlations in H3 exceed 1/sqrt(3), and those in
H4 exceed 1/sqrt(2), as follows directly from the inverses below.
Thus energy<=1/4 in H3, or energy<=11/36<1/3 in H4, forces ONE
actual opposite wall for the entire component. Every individual
half-tangent square is at most the component energy. Consequently

    h <= R^2 v_i^2 for every i,
    R^2=16/25 for H3,       R^2=1584/2209 for H4.       (1)

Here R^2=4e/(1+e)^2. These conclusions concern actual walls at the
same selected vertex; no coloring relaxation is being substituted.

## H3

In path order, with old bonds 5,3, the inverse unit Gram is

    C=[[3phi^2, 2phi^3, phi^3],
       [2phi^3, 4phi^2, 2phi^2],
       [phi^3, 2phi^2, phi+2]].

Assume energy<=1/4. Adjacent positive contacts with the common
opposite are impossible: deleting the remaining old node leaves a
nonpositive retained white triangle. If both endpoints have contacts
u,w>=1/2, put q for their own cost. Directly,

    q-(16/25)v_2^2
      =[(27phi+43)u^2+(9phi+34)w^2+(36phi+18)uw]/25
      >=(72phi+95)/100 > 1.

But h>=q-1 and (1) require this expression to be at most one.
Thus there is exactly one positive contact beta.

If that contact is at node1, already beta=1/2 gives

    q-(16/25)v_3^2=(199phi+263)/500 > 1.

At node2 the corresponding value, using direction3, is

    q-(16/25)v_3^2=(61phi+77)/125 > 1.

These expressions scale with beta^2 and are positive, so larger
contacts are impossible too. If the contact is at node3 and has
beta>=1/sqrt(2), using direction1 gives

    q-(16/25)v_1^2 >= (27phi+118)/150 > 1.

Only a native 3 contact at node3 remains. Its cost and smallest
directional derivative are

    q=(phi+2)/4<1,       v_1^2=phi^4/12.

Equation(1) gives h<=4phi^4/75<phi/4. If there is no spectator,
h=q-1<0. With exactly one positive spectator cost, Q>h since q<1,
contradicting the actual cross-component opposite obstruction.
With at least two, Q>=1/2, so h>=q-1/2=phi/4, again impossible.
The last strict radical comparison is 27phi>32. This proves H3>1/4
with arbitrary spectators, without enumerating them.

## H4

Use path bonds 5,3,3. The inverse unit Gram is

    C=[[8phi^4, 6phi^5, 4phi^5, 2phi^5],
       [6phi^5,12phi^4, 8phi^4, 4phi^4],
       [4phi^5, 8phi^4,16phi+12,8phi+6],
       [2phi^5, 4phi^4, 8phi+6,4phi^2]].

Assume energy<=11/36, so a common wall closes all four directions.
Any two contacts whose connecting old path omits an old node give
a nonpositive cycle at that node's retained finite endpoint. Thus
two contacts could only be at the two endpoints. But deleting old
node4 then leaves a four-node path with bond cosines at least
(1/2,phi/2,1/2), whose determinant is

    (5-4phi)/16<0.

It is not positive definite. Therefore only one contact is possible.

For a single contact at node j, test (1) on the direction i below.
At its smallest possible value beta=1/2, exact substitution gives:

|j|i|q-(1584/2209)v_i^2|
|---|---|---|
|1|4|(4938phi+3688)/2209|
|2|4|(7209phi+5334)/2209|
|3|1|(2500phi+2667)/2209|
|4|1|(625phi+1219)/2209|

Every entry is greater than one. The last and tightest comparison is
phi>198/125; this follows already from phi>8/5. Larger beta scales
each positive expression upward. Spectators can only increase h,
whereas h>=q-1 and (1) would require an entry at most one.
This proves the universal bound H4>11/36.

The proof imports the reviewed actual-wall alternatives,
orthogonal-opposite theorem and finite endpoint property. It has no
arithmetic hypothesis, external polytope classification or script
assumption. It leaves larger finite components and the global dimension
problem open. No additional manuscript or formalization is proposed
unless these local estimates lead to a material uniform bound.
