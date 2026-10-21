from __future__ import annotations

import json
import unittest

from scripts import certify_triangular_looijenga_h4_compiler as compiler


class TriangularLooijengaH4CompilerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.fresh = compiler.build_certificate()
        cls.archived = json.loads(compiler.OUTPUT.read_text(encoding="utf-8"))

    def test_archive_is_exact_and_all_checks_pass(self) -> None:
        self.assertEqual(self.archived, self.fresh)
        self.assertTrue(self.fresh["all_checks_pass"])
        self.assertFalse([name for name, passed in self.fresh["checks"].items() if not passed])

    def test_surface_free_triangle_and_branch_census_is_total(self) -> None:
        census = self.fresh["surface_free_census"]
        self.assertEqual(census["triangle_count"], 165)
        self.assertEqual(
            [row["count"] for row in census["triangle_sequence_profile"]],
            [60, 60, 20, 15, 10],
        )
        self.assertTrue(
            self.fresh["checks"][
                "all_triangles_have_plane_or_infinitely_near_blowup_template"
            ]
        )
        branches = self.fresh["branch_census"]
        self.assertEqual(branches["all_branch_count"], 885)
        self.assertEqual(branches["excluded_over_eight_facets"], 25)
        self.assertEqual(branches["surviving_decoration_branches"], 860)
        self.assertEqual(branches["surviving_marked_wall_sets"], 260)

    def test_four_positive_gram_types_have_exact_finite_volume_rays(self) -> None:
        rows = self.fresh["surviving_gram_types"]
        self.assertEqual([row["type_id"] for row in rows], ["T6a", "T6b", "T7", "T8"])
        self.assertEqual([row["facet_count"] for row in rows], [6, 6, 7, 8])
        self.assertEqual([row["marked_wall_set_count"] for row in rows], [120, 60, 60, 20])
        self.assertEqual(
            [
                (
                    row["cone_certificate"]["finite_vertex_count"],
                    row["cone_certificate"]["ideal_vertex_count"],
                    row["cone_certificate"]["hyperideal_ray_count"],
                )
                for row in rows
            ],
            [(5, 2, 0), (4, 3, 0), (5, 3, 0), (5, 4, 0)],
        )
        self.assertTrue(
            all(
                row["cone_certificate"]["decision"] == "finite_volume_noncompact"
                and all(row["checks"].values())
                for row in rows
            )
        )

    def test_inverse_recognizer_accepts_exactly_catalogued_representatives(self) -> None:
        rows = self.fresh["surviving_gram_types"]
        for row in rows:
            decision = compiler.recognize_raw(row["representative_walls"], rows)
            self.assertTrue(decision["carrier"])
            self.assertEqual(decision["type_id"], row["type_id"])
            self.assertEqual(decision["decision"], "finite_volume_noncompact")

        # Removing one T8 exceptional wall retains a raw vector family but is
        # not silently promoted to a complete triangular full-nef carrier.
        hostile = rows[-1]["representative_walls"][:-1]
        rejection = compiler.recognize_raw(hostile, rows)
        self.assertFalse(rejection["carrier"])
        self.assertIn(
            rejection["reason"],
            {
                "walls_do_not_span_H4_lattice",
                "unlisted_complete_nef_gram",
            },
        )


if __name__ == "__main__":
    unittest.main()
