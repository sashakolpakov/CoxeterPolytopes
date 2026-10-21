# Exact edge frontier of the 835-wall H22 branch

Every outgoing edge from the 23 previously certified finite vertices has now been classified against the existing 835 walls. These are **506 directed edge directions**, not an exhaustive list of all vertices of the chamber.

| Result | Directed directions |
|---|---:|
| Certified unclosed hyperbolic ray | 497 |
| Another existing finite endpoint | 7 |
| Actual ideal endpoint, with an unbounded Euclidean link | 2 |

The [frontier certificate](fresh_h22_pi7_edgewalk_20261013.json) verifies the first two rows. The [ideal-endpoint certificate](fresh_h22_pi7_edgewalk_20261013_ideal.json) resolves the two separately flagged cases. Thus neither ideal direction is mislabeled an open ray or a bounded cusp. The two unbounded links are subsequently repaired in a [coherent 837-wall extension](fresh_h22_pi7_cusp_pair_20261013.md); the classifications here refer to the frozen 835-wall chamber. The original [835-wall geometry](fresh_h22_incidence_prism_20261013_geometry.json) supplies the 23 initial vertices and the already known infinite-volume certificate.

For an incident positive Gram matrix A, let \(b_i=(h,w_i)\), where h is the strict timelike center. Define

\[
p=h-\sum_i(A^{-1}b)_i w_i,\qquad
u_j=\sum_i(A^{-1})_{ij}w_i.
\]

Then \(p\perp\nu_j\), \((\nu_j,w_i)=\delta_{ij}\), and \(\nu_j^2=(A^{-1})_{jj}>0\). The outgoing edge releasing wall j is \(p-s\nu_j\). For each of the 497 unclosed directions, the checker constructs

\[
x=p-(1+\epsilon)\sqrt{-p^2/\nu_j^2}\,\nu_j.
\]

All 21 retained wall pairings vanish identically. The other **814** pairings are strictly negative, \((x,h)<0\), and
\(x^2=(-p^2)(2\epsilon+\epsilon^2)>0\). Therefore the segment continues beyond the light cone without meeting any other wall. This is an exact existing-chamber obstruction, not merely failure to find an endpoint. In 476 directions \(\epsilon=1/16\) suffices; the other 21 use \(1/256\).

The seven finite directions have an exact positive edge parameter, negative endpoint norm, all nonincident inequalities strict, and a rigorously positive 22-wall Gram. Six directions form three pairs connecting vertices already in the list. The seventh reaches a newly identified rational finite vertex from saved vertex 22 by releasing wall 834 and meeting wall 38. Its 22 incident indices and full positive-Gram certificate are in that row. Its remaining outgoing directions are not processed here.

The two ideal directions release wall 21 at saved vertex 3, and wall 46 at vertex 4. Their primitive integer future null rays are archived. All 835 inequalities are checked. Each endpoint has 22 rational active walls, whose Gram has inertia \((20,0,2)\) and whose ambient span has dimension 21, plus one orthogonal positive-norm cap. The rational walls split into two affine components of ranks 8 and 12, each with an explicitly positive kernel vector. The cap supplies a one-sided additional Euclidean direction: its negative is a nonzero recession vector of the full link. The link is therefore unbounded, despite this being a genuine ideal endpoint.

A concrete constructive target is frontier row 0: saved vertex 0, release wall 1, retain cap 813 (new-cap index 1, candidate 140) and the listed 20 rational walls. This direction is among the exactly certified unclosed rays and is available for a new coupled endpoint cap. The certificates specify targets; they do not assert that suitable Coxeter caps exist.

All arithmetic inequalities use outward-rounded integer intervals with denominator \(10^{70}\). Exact projection identities handle vanishing pairings. The verified 835-wall root/angle data are inherited dependencies; this task does not independently rederive every old normal. A numerical [scout](fresh_h22_pi7_edgewalk_20261013_scout.json) only selected branches; every reported classification is subsequently certified.

Reproduce from the repository root:

```sh
PYTHONPATH=src:. python scripts/fresh_h22_pi7_edgewalk_20261013.py --verify
PYTHONPATH=src:. python scripts/fresh_h22_pi7_edgewalk_20261013_ideal.py --verify
```

The [independent source and formula review](fresh_h22_pi7_edgewalk_20261013_review.md) also passes both saved-result replays.

This is a finite frontier audit of a specified real Coxeter chamber. It is not a complete arithmetic edge-walking algorithm, a proof that all possible future caps fail, or a finite-volume result.
