# A blocker-decomposition screen for the 23 baseline H9/14 timeouts

## Screening theorem

Let \(H\) be a missing-edge graph on \(V=[14]\), let
\(C_1,\ldots,C_m\) be all inclusion-minimal vertex covers of \(H\), and let
\(\mathcal Q\) be the additional inclusion-minimal nonfaces of orders four
and five.  For nonempty \(J\subseteq[m]\), put

\[
U_J=\bigcup_{j\in J}C_j,
\qquad
\mathcal Q_J=\{Q\in\mathcal Q:Q\cap U_J=\varnothing\}.
\]

Write \(\tau_k(\mathcal Q_J;V\setminus U_J)\) for the number of \(k\)-sets
in \(V\setminus U_J\) meeting every member of \(\mathcal Q_J\).  Terms with
\(k<0\) are zero.  Then the number \(N_5\) of five-element transversals of
the complete nonface clutter is

\[
\boxed{
N_5=
\sum_{\varnothing\ne J\subseteq[m]\atop |U_J|\le5}
(-1)^{|J|+1}
\tau_{5-|U_J|}(\mathcal Q_J;V\setminus U_J).
}
\]

If the clutter has no transversal of order at most four, every five-element
transversal is inclusion-minimal.  Hence a pure eight-dimensional completion
must have no small transversal and must make the displayed integer equal the
scheduled facet number \(f_8\).

### Proof

Every transversal of the graph edges contains at least one minimal graph
cover \(C_i\).  Let \(X_i\) be the five-transversals containing \(C_i\).
An element of \(\bigcap_{j\in J}X_j\) contains \(U_J\).  Writing it uniquely
as \(U_J\sqcup S\), the remaining set \(S\) has order \(5-|U_J|\).  A member
of \(\mathcal Q\) meeting \(U_J\) is already hit, while every member disjoint
from \(U_J\) must be hit by \(S\).  Thus the cardinality of this intersection
is exactly the indicated \(\tau\)-term.  Inclusion--exclusion proves the
formula.  Finally, a nonminimal five-transversal contains a transversal of
order at most four; the converse implication is immediate.

Record 23 is the especially clean two-cover case in which this general formula
reduces to \(\tau_3(G)+\tau_2(Z)\).  The theorem shows that the same method is
available whenever the missing graph has only a few minimal covers, even when
their unions overlap and inclusion--exclusion is required.

## Exact screen of the remaining timeout shadows

The script `scripts/screen_h9_14_timeout_blockers.py` reconstructs the 23
timeouts in the frozen post-record-23 baseline, enumerates every
minimal graph cover over all \(2^{14}\) masks, and applies the strict inequality

\[
\sum_{v\in V}g_4(\operatorname{lk}v)\ge14
\]

to the full Dehn--Sommerville schedule.  In the table, `2:3,3:2` means three
minimal covers of order two and two of order three.  A term profile such as
`+4t2,-3t1,+2t0,-2t0` records the number, sign, and order of the exact
inclusion--exclusion terms; distinct terms of the same order are not presumed
equal.  The final column gives

\[
b(a)=\#\{\text{all missing independent five-sets}\},\quad
f_8(a),\quad
\sum_vg_4(\operatorname{lk}v).
\]

After missing-four shadows are removed, the number of additional minimal
five-generators is \(b(a)\) minus the shadow size.

| record | active graph + isolates | cover sizes | IE term profile | surviving \(a\) | \((b,f_8,\sum g_4)\) |
|---:|---|---|---|---:|---|
| 0 | 3v/2e + 11I | 1:1, 2:1 | +1t4,+1t3,-1t2 | 0–16 | (26+6a,180-2a,114-6a) |
| 1 | 4v/3e + 10I | 1:1, 3:1 | +1t4,+1t2,-1t1 | 0–13 | (22+6a,160-2a,96-6a) |
| 3 | 2K2 + 10I | 2:4 | +4t3,-4t2,+4t1,-3t1 | 0–13 | (20+6a,170-2a,96-6a) |
| 4 | 4v/3e + 10I | 2:3 | +3t3,-2t2,+1t1,-1t1 | 0–10 | (16+6a,150-2a,78-6a) |
| 5 | K2 + 3v/2e + 9I | 2:2, 3:2 | +2t3,+2t2,-1t2,-3t1,+4t0,-3t0 | 0–8 | (13+6a,142-2a,66-6a) |
| 6 | K3 + 11I | 2:3 | +3t3,+1t2,-3t2 | 0–10 | (16+6a,150-2a,78-6a) |
| 7 | 5v/4e + 9I | 2:2, 3:1 | +2t3,+1t2,-1t2,-1t1,+1t0,-1t0 | 0–7 | (12+6a,130-2a,60-6a) |
| 9 | 4v/4e + 10I | 2:2, 3:1 | +2t3,+1t2,-1t2,+1t1,-2t1 | 0–7 | (12+6a,130-2a,60-6a) |
| 12 | 4v/4e + 10I | 2:2 | +2t3,-1t1 | 0–7 | (12+6a,130-2a,60-6a) |
| 14 | 5v/5e + 9I | 2:1, 3:2 | +1t3,+2t2,-2t1,+1t0,-1t0 | 0–4 | (8+6a,110-2a,42-6a) |
| 16 | 4v/5e + 10I | 2:1, 3:2 | +1t3,+2t2,+1t1,-3t1 | 0–4 | (8+6a,110-2a,42-6a) |
| 17 | 5v/5e + 9I | 2:1, 3:2 | +1t3,+2t2,-2t1,+1t0,-1t0 | 0–4 | (8+6a,110-2a,42-6a) |
| 18 | 6v/5e + 8I | 2:1, 3:2 | +1t3,+2t2,-2t1 | 0–4 | (8+6a,110-2a,42-6a) |
| 31 | 5v/4e + 9I | 2:1, 3:3 | +1t3,+3t2,-4t1,+4t0,-3t0 | 0–5 | (9+6a,122-2a,48-6a) |
| 32 | (4v/3e)+K2+8I | 3:6 | +6t2,-7t1,+10t0,-8t0 | 0–3 | (6+6a,114-2a,36-6a) |
| 34 | K3+K2+9I | 3:6 | +6t2,+2t1,-9t1,+24t0,-22t0 | 0–3 | (6+6a,114-2a,36-6a) |
| 35 | 6v/5e + 8I | 3:5 | +5t2,-5t1,+6t0,-5t0 | 0–2 | (5+6a,102-2a,30-6a) |
| 37 | 5v/5e + 9I | 3:5 | +5t2,+1t1,-6t1,+10t0,-9t0 | 0–2 | (5+6a,102-2a,30-6a) |
| 43 | 6v/6e + 8I | 3:4 | +4t2,-3t1,+2t0,-2t0 | 0–1 | (4+6a,90-2a,24-6a) |
| 44 | 5v/6e + 9I | 3:4 | +4t2,-3t1,+4t0,-4t0 | 0–1 | (4+6a,90-2a,24-6a) |
| 46 | 5v/6e + 9I | 3:4 | +4t2,+1t1,-4t1,+3t0,-3t0 | 0–1 | (4+6a,90-2a,24-6a) |
| 47 | 6v/6e + 8I | 3:4 | +4t2,-3t1,+3t0,-3t0 | 0–1 | (4+6a,90-2a,24-6a) |
| 74 | K4 + 10I | 3:4 | +4t2,+4t1,-7t1 | 0–1 | (4+6a,90-2a,24-6a) |

All 23 shadows therefore have at most six minimal graph covers, so every
five-facet count decomposes into \(\tau_k\) terms with \(k\le4\).  The last
five records are the sharp first targets: strict link positivity leaves only
two schedules, all graph covers have order three, and only \(\tau_2\),
\(\tau_1\), and Boolean \(\tau_0\) terms occur.

## Record 43 outcome

The compressed `a=0` search for record 43 has 121 admissible active-row
patterns, 33 isolate-column signatures, and 1,073 replayed incidence states.
Its exact five-transversal support is

```text
58, 60–77 except 78, 79–83, 86, 90
```

and two compressed states attain the scheduled value 90.  One gives the four
missing five-faces displayed in `results/h9_14_record43_selection.md`.
The complete blocker, face table, all link audits, and a uniform rank-four
Gale sign certificate are frozen there. Thus, at this historical screening
stage, record 43 moved to the positive abstract-completion side rather than
producing a rejection.

Subsequent exact work resolves all five sharp targets from the bottom of the
table: records 43, 44, and 47 have explicit multiwedge realizations as
ordinary convex polytopes, while records 46 and 74 are independently
impossible. A later complete five-row quotient census also excludes record
37, and a later rational multiwedge realizes record 35. The current
record-34 certificate supplies one further exact exclusion. The current
theorem-facing partition is therefore

\[
121\text{ independently impossible}
\; +\;6\text{ explicit ordinary-polytopal selections}
\; +\;15\text{ theorem-undecided}.
\]

Only Coxeter labeling and Lorentzian Gram feasibility remain open for the
three Delzant multiwedge selections; they are no longer in the undecided
class. Record 35 is separately excluded as a compact Coxeter nerve.

## Reproduction

```bash
PYTHONPATH=src python3 scripts/screen_h9_14_timeout_blockers.py
PYTHONPATH=src python3 scripts/explore_h9_14_record43.py
```

The first command uses only exact graph enumeration and integer
Dehn--Sommerville/link arithmetic.  The second is a symmetry-compressed
discovery census; the final record-43 theorem does not trust its quotienting,
because the published explicit incidence structure is replayed directly over
all \(2^{14}\) masks by the independent certificate.
