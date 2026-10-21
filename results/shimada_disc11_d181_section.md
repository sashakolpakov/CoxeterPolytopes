# The exact physical \(D_{181}\) eigen-nef section

**Status.** Complete reduction of all primitive walls of the physical chamber
reached from \(D_{68}\) across its internal wall \((14,0)\), followed one
internal star farther. This is a local induced-chamber theorem, not a global
eigen-nef completeness, finite-volume, or compactness claim.

## 1. The physical state

The fixed object is the transported golden eigenspace of
`results/shimada_disc11_conjugate_marking.md`. Shimada's transition

\[
D_{68}\xrightarrow{(14,0)}D_{181}
\]

is applied to the chamber, not to the eigenspace. Thus every wall below is a
wall of the resulting physical copy of \(D_{181}\), restricted to the same
ten-dimensional real eigenspace.

An exact strict interior point is

\[
\begin{aligned}
p={}&\left(-\frac{359}{15}+\frac{31}{5}\phi,
 \frac{115}{6}-\frac{44}{3}\phi,
 3-3\phi,3-3\phi,5-5\phi,7-7\phi,\right.\\
&\left.11-11\phi,12-12\phi,15-15\phi,16-16\phi\right).
\end{aligned}
\]

Its trace square is

\[
q(p,p)=-\frac{11912}{25}+\frac{19643}{50}\phi>0,
\]

and all 26 primitive \(D_{181}\) wall functionals are strictly positive on
\(p\). Their exact minimum is

\[
\frac{481}{30}-\frac{99}{10}\phi>0.
\]

The external source is Shimada's 47 MB `disc11.txt`, with SHA-256

```text
fa39e5b59f4a4b2b1b7fe7632e5737d1734ee3eb7e7023ae26bf2b83387afb68
```

and 1,098 representative chambers.

## 2. Positive-proportional wall classes

Write \(v_{i,j}\) for member \(j\) of wall orbit \(i\), after physical
transport and restriction. Exact proportionality reduces the 26 primitive
walls to 17 distinct oriented half-spaces. The nine nonsingleton classes are

\[
\begin{array}{c|c}
\text{class representative}&\text{other member}\\ \hline
v_{0,0}&v_{0,1}=(\phi-1)v_{0,0}\\
v_{1,0}&v_{1,1}=\phi v_{1,0}\\
v_{2,0}&v_{6,0}=(\phi-1)v_{2,0}\\
v_{3,0}&v_{3,1}=(\phi-1)v_{3,0}\\
v_{4,0}&v_{4,1}=(\phi-1)v_{4,0}\\
v_{5,0}&v_{5,1}=\phi v_{5,0}\\
v_{7,0}&v_{7,1}=\phi v_{7,0}\\
v_{8,0}&v_{8,1}=\phi v_{8,0}\\
v_{9,0}&v_{9,1}=(\phi-1)v_{9,0}.
\end{array}
\]

All factors are positive at the distinguished embedding. The eight internal
walls

\[
v_{10,0},v_{10,1},v_{11,0},v_{12,0},v_{12,1},v_{13,0},v_{14,0},v_{15,0}
\]

are singleton classes.

The first ten Shimada wall orbits have `_minustwo` adjacency. Their 18
ambient K3 walls therefore restrict in nine pairs, exactly as predicted by
the universal even-golden two-root lift. One pair joins different Shimada
orbits, \(\{v_{2,0},v_{6,0}\}\); comparing only members inside one atlas
orbit would miss it.

## 3. Exact irredundant facets

The local section has precisely

\[
\boxed{11\text{ facets}=9\text{ root facets}+2\text{ internal facets}.}
\]

The root classes are

\[
\begin{gathered}
\{(0,0),(0,1)\},\ \{(1,0),(1,1)\},\ \{(2,0),(6,0)\},\\
\{(3,0),(3,1)\},\ \{(4,0),(4,1)\},\ \{(5,0),(5,1)\},\\
\{(7,0),(7,1)\},\ \{(8,0),(8,1)\},\ \{(9,0),(9,1)\}.
\end{gathered}
\]

The internal classes are

\[
\boxed{\{(10,0)\},\qquad\{(11,0)\}.}
\]

### Exact facet witnesses

No numerical reconstruction is used. If \(f\) is one of the eleven facet
functionals and \(H\) is the golden field Gram matrix, put

\[
x_f=p-\frac{f(p)}{fH^{-1}f^{\mathsf T}}H^{-1}f^{\mathsf T}. \tag{1}
\]

Every entry of \(x_f\) lies in \(\mathbf Q(\phi)\). Direct exact evaluation
proves \(f(x_f)=0\), strict positivity on every primitive wall outside the
target proportionality class, and \(q(x_f,x_f)>0\). Hence all eleven listed
half-spaces are genuine facets.

For the two internal facets, formula (1) gives

\[
\begin{aligned}
x_{10,0}={}&\frac1{180}\bigl(
-4789+1413\phi,2191-1862\phi,-2459+1313\phi,\\
&-2459+1313\phi,-2099+953\phi,-1739+593\phi,
-1019-127\phi,\\
&-839-307\phi,-299-847\phi,-119-1027\phi\bigr),
\end{aligned}
\]

with

\[
q(x_{10,0},x_{10,0})
=-\frac{226189}{1080}+\frac{13669}{60}\phi>0,
\]

and

\[
\begin{aligned}
x_{11,0}={}&\left(-\frac{80}{3}+8\phi,
\frac{52}{3}-\frac{40}{3}\phi,
3-3\phi,3-3\phi,5-5\phi,7-7\phi,\right.\\
&\left.11-11\phi,12-12\phi,15-15\phi,16-16\phi\right),
\end{aligned}
\]

with

\[
q(x_{11,0},x_{11,0})=-\frac{6988}{15}+\frac{1932}{5}\phi>0.
\]

The executable certificate emits all nine root-facet witnesses as exact
coordinate vectors too.

### Exact redundancy identities

Put

\[
f_i=v_{i,0}\quad(i=0,1,2,3,4,5,7,8,9),
\qquad a=v_{10,0},\qquad b=v_{11,0}.
\]

The six remaining half-spaces are redundant by the following identities. All
displayed coefficients are strictly positive at the identity embedding:

\[
v_{10,1}=(1+\phi)f_2+(\phi-1)a+2\phi b, \tag{2}
\]

\[
\begin{aligned}
v_{12,0}={}&(1+\phi)f_0+
\left(\frac12+\frac32\phi\right)f_1+
\left(\frac52+\frac92\phi\right)f_2+
\left(1-\frac12\phi\right)f_4\\
&+\left(\frac32+\frac52\phi\right)f_7+
(2+3\phi)f_8+f_9+
\left(1+\frac52\phi\right)a+(4+5\phi)b,
\end{aligned} \tag{3}
\]

\[
\begin{aligned}
v_{12,1}={}&\left(\frac13+\frac43\phi\right)f_0+
\left(\frac23+\frac43\phi\right)f_1+
\left(\frac{10}{3}+\frac{16}{3}\phi\right)f_2+
\left(\frac43-\frac23\phi\right)f_3\\
&+\left(-\frac13+\frac13\phi\right)f_5+
(2+2\phi)f_7+
\left(\frac83+\frac73\phi\right)f_8+f_9\\
&+\left(-\frac23+\frac{10}{3}\phi\right)a+
\left(4+\frac{20}{3}\phi\right)b,
\end{aligned} \tag{4}
\]

\[
\begin{aligned}
v_{13,0}={}&2\phi f_0+(2+2\phi)f_1+(1+\phi)f_2
+\phi f_3+\phi f_4\\
&+(2+2\phi)f_7+(2+2\phi)f_8+2\phi f_9+2\phi a+2\phi b,
\end{aligned} \tag{5}
\]

\[
\begin{aligned}
v_{14,0}={}&\frac53\phi f_0+
\left(\frac43+\frac43\phi\right)f_1+
\left(2+\frac{11}{3}\phi\right)f_2+
\frac23\phi f_3\\
&+\left(\frac13+\frac13\phi\right)f_5+
(2+2\phi)f_7+
\left(\frac73+\frac73\phi\right)f_8+
\phi f_9\\
&+\frac83\phi a+
\left(\frac73+5\phi\right)b,
\end{aligned} \tag{6}
\]

and

\[
v_{15,0}=f_2+\phi b. \tag{7}
\]

Equations (2)--(7) are checked coordinate by coordinate in
\(\mathbf Q(\phi)^{10}\). Together with the eleven strict face witnesses,
they prove the irredundant list without a floating-point facet decision.

## 4. One internal star farther

Both internal facets have been crossed using the exact opposite-side points

\[
z_f=x_f-\frac12(p-x_f).
\]

1. Facet \((10,0)\) has Shimada tag `_backto` and enters a physical copy of
   representative \(D_{30}\). The shared face has exactly one zero and no
   negative primitive-wall pairing in either chamber. Moreover
   \[
   q(z_{10,0},z_{10,0})
   =-\frac{1988609}{7200}+\frac{107631}{400}\phi>0,
   \]
   and \(z_{10,0}\) is strictly positive on all 34 primitive walls of that
   physical \(D_{30}\).

2. Facet \((11,0)\) has tag `_isom` and returns exactly to the original
   physical \(D_{68}\). Its face is the \(D_{68}\) wall \((14,0)\) crossed in
   the preceding step. Here
   \[
   q(z_{11,0},z_{11,0})
   =-\frac{11713}{25}+\frac{77603}{200}\phi>0,
   \]
   and \(z_{11,0}\) is strictly positive on all 28 primitive walls of
   \(D_{68}\).

Thus the internal adjacency of this local section is

\[
D_{30}^{\mathrm{phys}}
\xleftrightarrow{\ (10,0)\ }
D_{181}^{\mathrm{phys}}
\xleftrightarrow{\ (11,0)\ }
D_{68}.
\]

The six other internal primitive walls of representative \(D_{181}\) do not
cut the fixed eigensection and must not be traversed as section facets.

## 5. Scope and reproduction

The union of these three physical chambers is only a certified local patch.
In particular, the nine root classes are local facets rather than a global
simple-root list, and no volume or compactness conclusion follows.
The physical \(D_{30}\) has since been completely reduced and proved to be
the identical marked state reached through the \(D_{86}\) branch; see
`results/shimada_disc11_d30.md`.  Its new forward targets \(D_{16}\) and
\(D_{61}\) still require their own reductions.

Reproduce the result with

```sh
PYTHONPATH=src:. python3 scripts/analyze_shimada_disc11_d181.py \
  /path/to/disc11.txt
```

The script verifies the source checksum, all 17 proportionality classes, all
11 timelike face witnesses, all six positive dependence identities, and both
internal crossings.
