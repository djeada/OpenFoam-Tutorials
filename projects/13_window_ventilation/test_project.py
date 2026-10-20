#!/usr/bin/env python3
"""Regression tests for measurement definitions and rejection of bad output.

Run: python3 -m unittest discover -s . -p 'test_*.py' -v
Fixtures are synthetic parser inputs, never simulation results.
"""
import json
import contextlib
import io
import math
from pathlib import Path
import tempfile
import unittest

from analyze import history, read_table
from prepare_cases import geometry
from verify import scalar_field, verify_case


class MeasurementTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.case = Path(self.temp.name) / "wide_open"
        self.case.mkdir()
        (self.case / "case_info.json").write_text(json.dumps(dict(
            room_volume=1, room_temperature=293.15, end_time=2, write_interval=1,
            outdoor_temperature=273.15, radiator_temperature=313.15)))
        self.stats = self.case / "postProcessing/roomStats/0/volFieldValue.dat"
        self.flux = self.case / "postProcessing/windowFlux/0/surfaceFieldValue.dat"
        for path in (self.stats, self.flux):
            path.parent.mkdir(parents=True)
        self.stats.write_text("# Volume : 1\n# Time volAverage(c) volAverage(T)\n1 0.75 288.15\n2 0.5 283.15\n")
        self.flux.write_text("# Area : 0.12\n1 0.002\n2 0.004\n")

    def test_fresh_air_is_tracer_not_temperature(self):
        _, volume, rows = history(self.case)
        self.assertEqual(volume, 1)
        self.assertAlmostEqual(rows[-1]["fresh_air_percent"], 50)
        self.assertAlmostEqual(rows[-1]["mean_temperature_C"], 10)
        self.assertAlmostEqual(rows[-1]["exchange_ACH"], 7.2)
        self.assertAlmostEqual(rows[-1]["air_sensible_energy_decrease_J"], 12060)

    def test_inconsistent_time_histories_are_rejected(self):
        self.flux.write_text("1 0.002\n3 0.004\n")
        with self.assertRaisesRegex(ValueError, "unsynchronised"):
            history(self.case)

    def test_nonfinite_history_is_rejected(self):
        self.stats.write_text("1 nan 288.15\n")
        with self.assertRaisesRegex(ValueError, "invalid data"):
            read_table(self.case, "roomStats", "volFieldValue.dat", 3)

    def test_repeated_times_are_rejected(self):
        self.stats.write_text("1 0.9 290\n1 0.8 289\n")
        with self.assertRaisesRegex(ValueError, "non-increasing"):
            read_table(self.case, "roomStats", "volFieldValue.dat", 3)

    def completed_run(self):
        (self.case / "log.checkMesh").write_text("Mesh OK.\n")
        log = self.case / "log.buoyantBoussinesqPimpleFoam"
        log.write_text("Time = 2\nCourant Number mean: 0.1 max: 0.5\nEnd\n")
        final = self.case / "2"
        final.mkdir()
        for name, value in (("c", 0.5), ("T", 283.15), ("p_rgh", 0)):
            (final / name).write_text(f"internalField uniform {value};\n")
        return log, final

    def checked_case(self):
        with contextlib.redirect_stdout(io.StringIO()):
            return verify_case(self.case)

    def test_complete_output_passes_but_early_stop_fails(self):
        log, _ = self.completed_run()
        self.assertTrue(self.checked_case())
        log.write_text("Time = 1\nCourant Number mean: 0.1 max: 0.5\nEnd\n")
        self.assertFalse(self.checked_case())

    def test_local_tracer_overshoot_fails_even_with_valid_average(self):
        _, final = self.completed_run()
        path = final / "c"
        path.write_text("internalField nonuniform List<scalar> 2 (0.5 1.2);\n")
        self.assertFalse(self.checked_case())
        path.write_text("internalField nonuniform List<scalar> 3 (0.5 1.2);\n")
        with self.assertRaisesRegex(ValueError, "wrong scalar count"):
            scalar_field(path)


class GeometryTests(unittest.TestCase):
    def test_pane_displaces_air_and_window_plane_is_embedded(self):
        open_geo, open_volume = geometry(False, 0.08)
        tilted_geo, tilted_volume = geometry(True, 0.08)
        expected_solid = 0.04 * 1.2 * math.cos(math.radians(10)) + 0.0008
        self.assertAlmostEqual(open_volume, 1)
        self.assertAlmostEqual(tilted_volume, 1 - 0.1 * expected_solid)
        self.assertIn("In Surface{1}", open_geo)
        self.assertIn("In Surface{1}", tilted_geo)
        for patch in ("walls", "radiator", "outdoor", "frontAndBack"):
            self.assertIn(f'Physical Surface("{patch}")', tilted_geo)


if __name__ == "__main__":
    unittest.main()
