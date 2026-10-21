#!/usr/bin/env python3
"""Hostile replay of the H9/14 record-6 incidence exclusion.

This file intentionally does not import the production record-6 module.  It
rebuilds the fixed-base link census, the two labelled partition orbits, the
Dehn--Sommerville branch list, and every relative alignment in the ten
branches directly from subsets of an eleven-element set.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path
import shutil
import subprocess


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data/enumerations/h9_14_record6_incidence_no_go.json"

N = 11
V = tuple(range(N))
EDGES = tuple(combinations(V, 2))
TRIPLES = tuple(combinations(V, 3))
TRIPLE_POSITION = {triple: index for index, triple in enumerate(TRIPLES)}

EXPECTED_GENG_SHA256 = (
    "162cd507b93e74307ed9e88052da7bfb6f68e6cfeae6ccd2254c949c84ce0214"
)
EXPECTED_PRODUCTION_CERTIFICATE_SHA256 = (
    "cc0c3186b16cafc3c3a9572ea15f6d24880ad7bb3b18da40cf8c635d064f78ff"
)
EXPECTED_LINK_SURVIVOR_SHA256 = (
    "2d15881f939c4251fcd7e36bca901ade2fb66529a0b3481e266a51025fd9cb43"
)
EXPECTED_SUMMARY_SHA256 = (
    "42685136e71e5964620fdb7cd484784345f9e8bc63c4934d08dcd52e2267635e"
)

# (number of link facets, missing link triangles, g_4).
LINK_SIGNATURES = (
    (45, 2, 1),
    (48, 1, 1),
    (49, 1, 2),
    (50, 1, 3),
    (52, 0, 2),
    (53, 0, 3),
    (54, 0, 4),
    (55, 0, 5),
)

EXPECTED_LINK_HISTOGRAM = (
    ((45, 2, 1), 2),
    ((48, 1, 1), 2),
    ((49, 1, 2), 10),
    ((50, 1, 3), 4),
    ((52, 0, 2), 6),
    ((53, 0, 3), 12),
    ((54, 0, 4), 5),
    ((55, 0, 5), 1),
)

# (three sorted link signatures, isolate-only missing A4 count, |Z|).
EXPECTED_NUMERICAL_BRANCHES = (
    (((45, 2, 1), (45, 2, 1), (45, 2, 1)), 0, 3),
    (((45, 2, 1), (45, 2, 1), (45, 2, 1)), 1, 1),
    (((45, 2, 1), (45, 2, 1), (48, 1, 1)), 0, 2),
    (((45, 2, 1), (45, 2, 1), (48, 1, 1)), 1, 0),
    (((45, 2, 1), (45, 2, 1), (49, 1, 2)), 0, 1),
    (((45, 2, 1), (45, 2, 1), (50, 1, 3)), 0, 0),
    (((45, 2, 1), (45, 2, 1), (52, 0, 2)), 0, 0),
    (((45, 2, 1), (48, 1, 1), (48, 1, 1)), 0, 1),
    (((45, 2, 1), (48, 1, 1), (49, 1, 2)), 0, 0),
    (((48, 1, 1), (48, 1, 1), (48, 1, 1)), 0, 0),
)

# name, facet-count candidates, link-compatible candidates, candidates with
# every isolate triple covered, candidates with the required a0, missing-
# triple histogram, and (after triple coverage) isolate-A4 histogram.
EXPECTED_BRANCH_SUMMARY = (
    ("EEE_z0_a0_0", 54, 54, 36, 0, ((0, 36), (1, 18)), ((1, 36),)),
    ("AHH_z0_a0_0", 0, 0, 0, 0, (), ()),
    ("EEH_z1_a0_0", 2_790, 0, 0, 0, (), ()),
    ("EFH_z0_a0_0", 0, 0, 0, 0, (), ()),
    ("EHH_z2_a0_0", 25_086, 0, 0, 0, (), ()),
    ("FHH_z1_a0_0", 1_953, 1_440, 0, 0, ((1, 1_440),), ()),
    ("GHH_z0_a0_0", 0, 0, 0, 0, (), ()),
    ("HHH_z3_a0_0", 780, 780, 0, 0, ((1, 780),), ()),
    ("EHH_z0_a0_1", 60, 60, 0, 0, ((1, 60),), ()),
    ("HHH_z1_a0_1", 1_130, 0, 0, 0, (), ()),
)


def mask(vertices) -> int:
    return sum(1 << vertex for vertex in vertices)


TRIPLE_VERTEX_MASKS = tuple(mask(triple) for triple in TRIPLES)
PAIR_MASKS = tuple(mask(edge) for edge in EDGES)
FOUR_MASKS = tuple(mask(four) for four in combinations(V, 4))
FOUR_TRIPLE_MASKS = tuple(
    sum(
        1 << index
        for index, triple_mask in enumerate(TRIPLE_VERTEX_MASKS)
        if triple_mask & four_mask == triple_mask
    )
    for four_mask in FOUR_MASKS
)
DISJOINT_TRIPLES = tuple(
    sum(
        1 << index
        for index, triple_mask in enumerate(TRIPLE_VERTEX_MASKS)
        if not subset & triple_mask
    )
    for subset in range(1 << N)
)
TRIPLE_NEIGHBOURS = tuple(
    sum(
        1 << other
        for other, other_mask in enumerate(TRIPLE_VERTEX_MASKS)
        if other != index and (triple_mask | other_mask).bit_count() == 4
    )
    for index, triple_mask in enumerate(TRIPLE_VERTEX_MASKS)
)
EDGE_COBoundaries = tuple(
    sum(
        1 << index
        for index, triple in enumerate(TRIPLES)
        if first in triple and second in triple
    )
    for first, second in EDGES
)


def decode_graph6(word: str) -> tuple[int, tuple[tuple[int, int], ...]]:
    """Decode the short graph6 form used by ``geng -q -t 10``."""

    codes = [ord(character) - 63 for character in word]
    if not codes or not 0 <= codes[0] <= 62 or any(not 0 <= code <= 63 for code in codes):
        raise ValueError("unsupported graph6 word")
    order = codes[0]
    bit_count = order * (order - 1) // 2
    if len(codes) != 1 + (bit_count + 5) // 6:
        raise ValueError("wrong graph6 payload length")
    bits = []
    for code in codes[1:]:
        bits.extend((code >> shift) & 1 for shift in range(5, -1, -1))
    if any(bits[bit_count:]):
        raise ValueError("nonzero graph6 padding")
    result = []
    position = 0
    for upper in range(1, order):
        for lower in range(upper):
            if bits[position]:
                result.append((lower, upper))
            position += 1
    return order, tuple(result)


def ridge_switching_graph(edges: tuple[tuple[int, int], ...]) -> bool:
    """Test triangle-free and induced-2K2-free without production helpers."""

    edge_set = set(edges)
    if any(
        all(tuple(sorted(pair)) in edge_set for pair in combinations(triple, 2))
        for triple in combinations(range(10), 3)
    ):
        return False
    for four in combinations(range(10), 4):
        induced = [
            tuple(sorted(pair))
            for pair in combinations(four, 2)
            if tuple(sorted(pair)) in edge_set
        ]
        if len(induced) == 2 and set(induced[0]).isdisjoint(induced[1]):
            return False
    return True


def switching_two_graph(edges: tuple[tuple[int, int], ...]) -> int:
    """Add a switch base 0 and form delta(X) on eleven vertices."""

    lifted = {(first + 1, second + 1) for first, second in edges}
    result = 0
    for index, triple in enumerate(TRIPLES):
        parity = sum(tuple(sorted(pair)) in lifted for pair in combinations(triple, 2))
        result |= (parity & 1) << index
    return result


def connected_facet_dual(two_graph: int) -> bool:
    remaining = two_graph
    if not remaining:
        return False
    first = (remaining & -remaining).bit_length() - 1
    reached = 1 << first
    frontier = [first]
    while frontier:
        current = frontier.pop()
        fresh = TRIPLE_NEIGHBOURS[current] & two_graph & ~reached
        while fresh:
            bit = fresh & -fresh
            reached |= bit
            frontier.append(bit.bit_length() - 1)
            fresh ^= bit
    return reached == two_graph


@lru_cache(maxsize=None)
def link_signature(two_graph: int) -> tuple[int, int, int] | None:
    """Recompute all local incidence necessities from the link facets."""

    if any((two_graph & four).bit_count() not in (0, 2) for four in FOUR_TRIPLE_MASKS):
        return None
    if any(not two_graph & DISJOINT_TRIPLES[pair] for pair in PAIR_MASKS):
        return None

    face = tuple(
        subset == 0 or bool(two_graph & DISJOINT_TRIPLES[subset])
        for subset in range(1 << N)
    )
    minimal_nonfaces = []
    for subset in range(1, 1 << N):
        if face[subset]:
            continue
        deletions_are_faces = True
        remaining = subset
        while remaining:
            bit = remaining & -remaining
            if not face[subset ^ bit]:
                deletions_are_faces = False
                break
            remaining ^= bit
        if not deletions_are_faces:
            continue
        if subset.bit_count() > 5:
            return None
        minimal_nonfaces.append(subset)
    if not minimal_nonfaces or not connected_facet_dual(two_graph):
        return None

    missing_triangles = sum(value.bit_count() == 3 for value in minimal_nonfaces)
    if missing_triangles not in (0, 1, 2, 3):
        return None
    four_faces = sum(face[value] for value in FOUR_MASKS)
    g4 = four_faces - 325 + 6 * missing_triangles
    # Macaulay raises of g3=4-q in degree three are 5,3,2,1.
    if not 1 <= g4 <= (5, 3, 2, 1)[missing_triangles]:
        return None
    return two_graph.bit_count(), missing_triangles, g4


def transversal_two_graph(blocks) -> int:
    result = 0
    first, second, third = (tuple(block) for block in blocks)
    for a in first:
        for b in second:
            for c in third:
                result |= 1 << TRIPLE_POSITION[tuple(sorted((a, b, c)))]
    return result


def e_orbit() -> tuple[int, ...]:
    """Enumerate 3+4+4 partitions by rooting one equal block at min(rest)."""

    result = []
    for small in combinations(V, 3):
        rest = tuple(vertex for vertex in V if vertex not in small)
        for companions in combinations(rest[1:], 3):
            first_large = (rest[0],) + companions
            second_large = tuple(vertex for vertex in rest if vertex not in first_large)
            result.append(transversal_two_graph((small, first_large, second_large)))
    if len(result) != 5_775 or len(set(result)) != 5_775:
        raise AssertionError("independent E-orbit enumeration changed")
    return tuple(result)


def h_orbit() -> tuple[int, ...]:
    """Enumerate 3+3+5 partitions by rooting one equal block at min(rest)."""

    result = []
    for large in combinations(V, 5):
        rest = tuple(vertex for vertex in V if vertex not in large)
        for companions in combinations(rest[1:], 2):
            first_small = (rest[0],) + companions
            second_small = tuple(vertex for vertex in rest if vertex not in first_small)
            result.append(transversal_two_graph((first_small, second_small, large)))
    if len(result) != 4_620 or len(set(result)) != 4_620:
        raise AssertionError("independent H-orbit enumeration changed")
    return tuple(result)


E0 = transversal_two_graph((range(3), range(3, 7), range(7, 11)))
H0 = transversal_two_graph((range(3), range(3, 6), range(6, 11)))


def edge_coboundary(edge_indices) -> int:
    result = 0
    for index in edge_indices:
        result ^= EDGE_COBoundaries[index]
    return result


def derive_numerical_branches():
    result = []
    for signatures in combinations_with_replacement(LINK_SIGNATURES, 3):
        link_facets = sum(signature[0] for signature in signatures)
        through_special = sum(signature[1] for signature in signatures)
        a0 = 0
        while True:
            z_size = 150 - link_facets - 2 * (through_special + a0)
            if z_size < 0:
                break
            result.append((tuple(signatures), a0, z_size))
            a0 += 1
    return tuple(sorted(result))


def isolate_nonfaces(ys, z_edge_indices, size: int) -> tuple[tuple[int, ...], ...]:
    selected_y = 0
    for two_graph in ys:
        selected_y |= two_graph
    z_masks = tuple(PAIR_MASKS[index] for index in z_edge_indices)
    result = []
    for chosen in combinations(V, size):
        chosen_mask = mask(chosen)
        covered_by_y = bool(selected_y & DISJOINT_TRIPLES[chosen_mask])
        covered_by_z = any(not chosen_mask & edge_mask for edge_mask in z_masks)
        if not covered_by_y and not covered_by_z:
            result.append(chosen)
    return tuple(result)


@dataclass(frozen=True)
class BranchResult:
    name: str
    link_types: str
    isolate_missing_four_count: int
    z_size: int
    target_signature: tuple[int, int, int]
    facet_count_candidates: int
    link_compatible_candidates: int
    independent_triple_compatible_candidates: int
    expected_a0_compatible_candidates: int
    missing_triple_histogram: tuple[tuple[int, int], ...]
    isolate_four_histogram: tuple[tuple[int, int], ...]

    def frozen_summary(self):
        return (
            self.name,
            self.facet_count_candidates,
            self.link_compatible_candidates,
            self.independent_triple_compatible_candidates,
            self.expected_a0_compatible_candidates,
            self.missing_triple_histogram,
            self.isolate_four_histogram,
        )

    def to_record(self) -> dict[str, object]:
        return {
            "name": self.name,
            "link_types": self.link_types,
            "isolate_missing_four_count": self.isolate_missing_four_count,
            "z_size": self.z_size,
            "target_signature": list(self.target_signature),
            "facet_count_candidates": self.facet_count_candidates,
            "link_compatible_candidates": self.link_compatible_candidates,
            "independent_triple_compatible_candidates": (
                self.independent_triple_compatible_candidates
            ),
            "expected_a0_compatible_candidates": self.expected_a0_compatible_candidates,
            "missing_triple_histogram": {
                str(key): value for key, value in self.missing_triple_histogram
            },
            "isolate_four_histogram": {
                str(key): value for key, value in self.isolate_four_histogram
            },
        }


def scan_branch(
    *,
    name: str,
    link_types: str,
    expected_a0: int,
    base: int,
    orbit: tuple[int, ...],
    z_supports: tuple[tuple[int, ...], ...],
    target_signature: tuple[int, int, int],
) -> BranchResult:
    facet_count_candidates = 0
    link_candidates = []
    for second in orbit:
        for z_edges in z_supports:
            third = base ^ second ^ edge_coboundary(z_edges)
            if third.bit_count() != target_signature[0]:
                continue
            facet_count_candidates += 1
            if link_signature(third) == target_signature:
                link_candidates.append((second, third, z_edges))

    triple_histogram = Counter()
    four_histogram = Counter()
    triple_compatible = 0
    expected_compatible = 0
    for second, third, z_edges in link_candidates:
        missing_triples = isolate_nonfaces((base, second, third), z_edges, 3)
        triple_histogram[len(missing_triples)] += 1
        if missing_triples:
            continue
        triple_compatible += 1
        missing_fours = isolate_nonfaces((base, second, third), z_edges, 4)
        four_histogram[len(missing_fours)] += 1
        expected_compatible += len(missing_fours) == expected_a0

    return BranchResult(
        name=name,
        link_types=link_types,
        isolate_missing_four_count=expected_a0,
        z_size=len(z_supports[0]),
        target_signature=target_signature,
        facet_count_candidates=facet_count_candidates,
        link_compatible_candidates=len(link_candidates),
        independent_triple_compatible_candidates=triple_compatible,
        expected_a0_compatible_candidates=expected_compatible,
        missing_triple_histogram=tuple(sorted(triple_histogram.items())),
        isolate_four_histogram=tuple(sorted(four_histogram.items())),
    )


def scan_hhh_z3(h_links: tuple[int, ...]) -> BranchResult:
    # A three-edge graph has a unique coboundary here: two such graphs could
    # differ by a nonzero cut only if their symmetric difference had at least
    # ten edges, whereas it has at most six.
    three_edge_coboundaries = {
        edge_coboundary(edges): edges for edges in combinations(range(55), 3)
    }
    if len(three_edge_coboundaries) != 26_235:
        raise AssertionError("three-edge coboundary map is not injective")

    unordered = []
    for first_index, first in enumerate(h_links):
        for second_index in range(first_index, len(h_links)):
            delta = H0 ^ first ^ h_links[second_index]
            z_edges = three_edge_coboundaries.get(delta)
            if z_edges is not None:
                unordered.append((first_index, second_index, z_edges))

    triple_histogram = Counter()
    four_histogram = Counter()
    ordered_count = 0
    triple_compatible = 0
    expected_compatible = 0
    for first_index, second_index, z_edges in unordered:
        multiplicity = 1 if first_index == second_index else 2
        ordered_count += multiplicity
        ys = (H0, h_links[first_index], h_links[second_index])
        missing_triples = isolate_nonfaces(ys, z_edges, 3)
        triple_histogram[len(missing_triples)] += multiplicity
        if missing_triples:
            continue
        triple_compatible += multiplicity
        missing_fours = isolate_nonfaces(ys, z_edges, 4)
        four_histogram[len(missing_fours)] += multiplicity
        expected_compatible += multiplicity * (not missing_fours)

    return BranchResult(
        name="HHH_z3_a0_0",
        link_types="H,H,H",
        isolate_missing_four_count=0,
        z_size=3,
        target_signature=(45, 2, 1),
        facet_count_candidates=ordered_count,
        link_compatible_candidates=ordered_count,
        independent_triple_compatible_candidates=triple_compatible,
        expected_a0_compatible_candidates=expected_compatible,
        missing_triple_histogram=tuple(sorted(triple_histogram.items())),
        isolate_four_histogram=tuple(sorted(four_histogram.items())),
    )


def link_census(e_links: set[int], h_links: set[int]) -> dict[str, object]:
    executable = shutil.which("geng")
    if executable is None:
        raise RuntimeError("nauty geng is required for the hostile link replay")
    completed = subprocess.run(
        [executable, "-q", "-t", "10"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    transcript_hash = sha256(completed.stdout).hexdigest()
    if transcript_hash != EXPECTED_GENG_SHA256:
        raise AssertionError("triangle-free graph transcript changed")

    survivors = []
    signatures = []
    two_graphs = []
    ridge_count = 0
    for raw_word in completed.stdout.splitlines():
        word = raw_word.decode("ascii")
        order, edges = decode_graph6(word)
        if order != 10:
            raise AssertionError("mixed-order geng transcript")
        if not ridge_switching_graph(edges):
            continue
        ridge_count += 1
        two_graph = switching_two_graph(edges)
        signature = link_signature(two_graph)
        if signature is not None:
            survivors.append(word)
            signatures.append(signature)
            two_graphs.append(two_graph)

    survivor_hash = sha256(
        "".join(word + "\n" for word in survivors).encode()
    ).hexdigest()
    histogram = tuple(sorted(Counter(signatures).items()))
    if (
        len(completed.stdout.splitlines()) != 12_172
        or ridge_count != 308
        or len(survivors) != 42
        or histogram != EXPECTED_LINK_HISTOGRAM
    ):
        raise AssertionError("independent fixed-base link census changed")
    if (
        EXPECTED_LINK_SURVIVOR_SHA256 != "TO_BE_FROZEN"
        and survivor_hash != EXPECTED_LINK_SURVIVOR_SHA256
    ):
        raise AssertionError("fixed-base link survivor stream changed")

    e_survivors = [
        two_graph
        for two_graph, signature in zip(two_graphs, signatures)
        if signature == (48, 1, 1)
    ]
    h_survivors = [
        two_graph
        for two_graph, signature in zip(two_graphs, signatures)
        if signature == (45, 2, 1)
    ]
    if len(e_survivors) != 2 or not set(e_survivors) <= e_links:
        raise AssertionError("E source-link orbit is incomplete")
    if len(h_survivors) != 2 or not set(h_survivors) <= h_links:
        raise AssertionError("H source-link orbit is incomplete")
    return {
        "triangle_free_unlabelled_graphs": 12_172,
        "triangle_and_induced_2k2_free_graphs": ridge_count,
        "survivor_count": len(survivors),
        "survivor_stream_sha256": survivor_hash,
        "signature_histogram": [
            {"signature": list(signature), "multiplicity": multiplicity}
            for signature, multiplicity in histogram
        ],
        "E_survivors_in_full_orbit": len(e_survivors),
        "H_survivors_in_full_orbit": len(h_survivors),
    }


def ten_branch_replay(e_links: tuple[int, ...], h_links: tuple[int, ...]):
    z0 = ((),)
    z1 = tuple((index,) for index in range(55))
    z2 = tuple(combinations(range(55), 2))
    branches = (
        scan_branch(
            name="EEE_z0_a0_0", link_types="E,E,E", expected_a0=0,
            base=E0, orbit=e_links, z_supports=z0,
            target_signature=(48, 1, 1),
        ),
        scan_branch(
            name="AHH_z0_a0_0", link_types="H,H,A", expected_a0=0,
            base=H0, orbit=h_links, z_supports=z0,
            target_signature=(52, 0, 2),
        ),
        scan_branch(
            name="EEH_z1_a0_0", link_types="H,E,E", expected_a0=0,
            base=H0, orbit=e_links, z_supports=z1,
            target_signature=(48, 1, 1),
        ),
        scan_branch(
            name="EFH_z0_a0_0", link_types="H,E,F", expected_a0=0,
            base=H0, orbit=e_links, z_supports=z0,
            target_signature=(49, 1, 2),
        ),
        scan_branch(
            name="EHH_z2_a0_0", link_types="H,H,E", expected_a0=0,
            base=H0, orbit=h_links, z_supports=z2,
            target_signature=(48, 1, 1),
        ),
        scan_branch(
            name="FHH_z1_a0_0", link_types="H,H,F", expected_a0=0,
            base=H0, orbit=h_links, z_supports=z1,
            target_signature=(49, 1, 2),
        ),
        scan_branch(
            name="GHH_z0_a0_0", link_types="H,H,G", expected_a0=0,
            base=H0, orbit=h_links, z_supports=z0,
            target_signature=(50, 1, 3),
        ),
        scan_hhh_z3(h_links),
        scan_branch(
            name="EHH_z0_a0_1", link_types="H,H,E", expected_a0=1,
            base=H0, orbit=h_links, z_supports=z0,
            target_signature=(48, 1, 1),
        ),
        scan_branch(
            name="HHH_z1_a0_1", link_types="H,H,H", expected_a0=1,
            base=H0, orbit=h_links, z_supports=z1,
            target_signature=(45, 2, 1),
        ),
    )
    if tuple(branch.frozen_summary() for branch in branches) != EXPECTED_BRANCH_SUMMARY:
        raise AssertionError("independent ten-branch failure histogram changed")
    return branches


def canonical_digest(document: dict[str, object]) -> str:
    return sha256(
        json.dumps(document, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def manifest_branch_projection(row: dict[str, object]):
    return (
        row["name"],
        row["facet_count_candidates"],
        row["link_compatible_candidates"],
        row["independent_triple_compatible_candidates"],
        row["expected_a0_compatible_candidates"],
        tuple(sorted((int(key), value) for key, value in row["missing_triple_histogram"].items())),
        tuple(sorted((int(key), value) for key, value in row["isolate_four_histogram"].items())),
    )


def audit(manifest_path: Path = MANIFEST) -> dict[str, object]:
    numerical = derive_numerical_branches()
    if numerical != EXPECTED_NUMERICAL_BRANCHES:
        raise AssertionError("independent DS branch derivation changed")

    e_links = e_orbit()
    h_links = h_orbit()
    if link_signature(E0) != (48, 1, 1) or link_signature(H0) != (45, 2, 1):
        raise AssertionError("canonical E/H link signature changed")
    census = link_census(set(e_links), set(h_links))
    branches = ten_branch_replay(e_links, h_links)

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload = dict(manifest)
    claimed_certificate_hash = payload.pop("certificate_sha256")
    actual_certificate_hash = canonical_digest(payload)
    if (
        claimed_certificate_hash != actual_certificate_hash
        or actual_certificate_hash != EXPECTED_PRODUCTION_CERTIFICATE_SHA256
    ):
        raise AssertionError("production certificate digest changed")
    production_projection = tuple(
        manifest_branch_projection(row) for row in manifest["branch_audits"]
    )
    hostile_projection = tuple(branch.frozen_summary() for branch in branches)
    if production_projection != hostile_projection:
        raise AssertionError("hostile and production branch counts disagree")

    summary: dict[str, object] = {
        "kind": "h9_14_record6_independent_hostile_replay",
        "imports_production_module": False,
        "link_census": census,
        "source_link_orbits": {
            "E_3_plus_4_plus_4_labelled_count": len(e_links),
            "H_3_plus_3_plus_5_labelled_count": len(h_links),
        },
        "numerical_branches": [
            {
                "link_signatures": [list(signature) for signature in signatures],
                "isolate_missing_four_count": a0,
                "z_size": z_size,
            }
            for signatures, a0, z_size in numerical
        ],
        "branch_audits": [branch.to_record() for branch in branches],
        "production_certificate_sha256": actual_certificate_hash,
        "production_counts_match": True,
        "all_ten_branches_excluded": all(
            branch.expected_a0_compatible_candidates == 0 for branch in branches
        ),
    }
    digest = canonical_digest(summary)
    if EXPECTED_SUMMARY_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_SUMMARY_SHA256:
        raise AssertionError("hostile replay summary digest changed")
    summary["summary_sha256"] = digest
    return summary


def main() -> None:
    print(json.dumps(audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
