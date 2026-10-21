# Conditional exclusion of a 30-vertex terminal graph

This is an independent abstract graph review. It does not reconstruct the
geometric terminal graph or invoke the Andrásfai–Erdős–Sós theorem from a
newly checked source. The following four properties are assumptions:

1. The graph T has 30 vertices, is K5-free, and has minimum degree at least22.
2. Its vertex set partitions into five induced copies of K_{2,2,2}.
3. At most two vertices of T are adjacent to all six vertices of any block.
4. The AES theorem applies in its form: a K5-free graph on n vertices with
   minimum degree greater than8n/11 is four-partite.

These assumptions are incompatible. Indeed22>240/11, so take a proper
four-coloring. Each color class has size at most8, because a vertex has
at least22 neighbors outside its own class. The four class sizes, in
increasing order, are therefore (6,8,8,8) or (7,7,8,8).

A vertex in a class of size8 is adjacent to every vertex outside that
class. A vertex in a class of size7 misses at most one outside vertex,
and one in a class of size6 misses at most two.

Every K_{2,2,2} block uses at least three colors. If it uses three, its
three opposite pairs are monochromatic and it omits one class. An
omitted class of size8 supplies eight common neighbors by the preceding
universality. If the omitted class has size6, the three used classes
all have size8, and all six vertices of the omitted class are common
neighbors. If it has size7, the only block vertices which can miss a
vertex of that omitted class are the two in the other class of size7.
Each misses at most one, leaving at least five common neighbors. Every
case contradicts assumption3.

Thus every block uses all four colors. Colors used on distinct opposite
pairs are disjoint, since such vertices are adjacent. Exactly one
opposite pair is bichromatic. Its two vertices are nonadjacent, so
neither can belong to a class of size8. Consequently the class sizes
are (7,7,8,8), and this pair uses one vertex of each size7 class. Each
of the remaining two monochromatic pairs is adjacent to both members
of that pair; the pairs must therefore use the two size8 classes, one
pair in each. Each block uses color counts (1,1,2,2). Five disjoint
blocks would use counts (5,5,10,10), contradicting the actual class
sizes (7,7,8,8).

The arithmetic and coloring argument require no enumeration. Turning
this conditional exclusion into a geometric lower bound still requires
an intrinsic proof of assumptions1–3 and an explicit supported AES input.
No claim about the finite-volume Coxeter ambient57 frontier follows.
