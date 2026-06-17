#!/usr/bin/env python3

import unittest
import sys
import os
import io
import json
import gzip
import tempfile
import shutil
from contextlib import redirect_stderr
from types import SimpleNamespace
from unittest.mock import patch

# Add parent directory to path to import the module under test
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import meta_analysis
from meta_analysis import (
    ordered_groups,
    loo_studies,
    loo_header_cols,
    Study,
)


def study_ns(**attrs):
    """Lightweight stand-in for Study exposing only the attributes used by the
    grouping helpers (cohort / population / name)."""
    attrs.setdefault("cohort", None)
    attrs.setdefault("population", None)
    return SimpleNamespace(**attrs)


class TestGroupingHelpers(unittest.TestCase):
    """Unit tests for the pure grouping helpers."""

    def test_ordered_groups_all_tagged(self):
        studs = [study_ns(cohort="A"), study_ns(cohort="B"), study_ns(cohort="A")]
        self.assertEqual(ordered_groups(studs, "cohort"), ["A", "B"])

    def test_ordered_groups_none_tagged(self):
        studs = [study_ns(), study_ns(), study_ns()]
        self.assertEqual(ordered_groups(studs, "cohort"), [])

    def test_ordered_groups_partial(self):
        # Only non-missing values, deterministic config order, field-less skipped
        studs = [study_ns(population="EUR"), study_ns(), study_ns(population="AFR"),
                 study_ns(population="EUR")]
        self.assertEqual(ordered_groups(studs, "population"), ["EUR", "AFR"])

    def test_loo_studies_drops_group_keeps_none(self):
        s0 = study_ns(cohort="A")
        s1 = study_ns(cohort="A")
        s2 = study_ns(cohort="B")
        s3 = study_ns(cohort=None)  # field-less, always retained
        studs = [s0, s1, s2, s3]
        next_var = ["v0", "v1", "v2", "v3"]  # all present (non-None)
        kept = loo_studies(next_var, studs, "cohort", "A")
        kept_studies = [pair[0] for pair in kept]
        self.assertEqual(kept_studies, [s2, s3])  # A dropped, B and field-less kept

    def test_loo_studies_skips_absent_variant(self):
        s0 = study_ns(cohort="A")
        s1 = study_ns(cohort="B")
        studs = [s0, s1]
        next_var = [None, "v1"]  # s0 has no variant at this position
        kept = loo_studies(next_var, studs, "cohort", "B")
        self.assertEqual(kept, [])  # B dropped, s0 absent -> nothing

    def test_loo_header_cols_no_het(self):
        cols = loo_header_cols("cohort_A", ["n"], is_het_test=False)
        self.assertEqual(cols, [
            "leave_cohort_A_N",
            "leave_cohort_A_n_meta_beta",
            "leave_cohort_A_n_meta_sebeta",
            "leave_cohort_A_n_meta_p",
            "leave_cohort_A_n_meta_mlogp",
        ])

    def test_loo_header_cols_with_het(self):
        cols = loo_header_cols("FINNGEN", ["inv_var"], is_het_test=True)
        self.assertIn("leave_FINNGEN_inv_var_meta_het_p", cols)
        self.assertEqual(len(cols), 1 + 5)  # N + 5 cols for one method with het


class TestStudyMetadataParsing(unittest.TestCase):
    """The cohort / population fields are optional config metadata, not data
    columns, so they must parse without touching the data file header."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def _sumstats(self, name):
        path = os.path.join(self.temp_dir, name)
        with gzip.open(path, "wt") as f:
            f.write("CHR\tPOS\tREF\tALT\tbeta\tpval\n")
            f.write("1\t1000\tA\tG\t0.1\t0.01\n")
        return path

    def _base_conf(self, fpath):
        return {
            "name": "S1", "file": fpath, "n_cases": 100, "n_controls": 200,
            "chr": "CHR", "pos": "POS", "ref": "REF", "alt": "ALT",
            "effect": "beta", "effect_type": "beta", "pval": "pval",
        }

    def test_metadata_present(self):
        conf = self._base_conf(self._sumstats("s1.gz"))
        conf["cohort"] = "A"
        conf["population"] = "EUR"
        s = Study(conf)
        self.assertEqual(s.cohort, "A")
        self.assertEqual(s.population, "EUR")

    def test_metadata_absent_is_none(self):
        s = Study(self._base_conf(self._sumstats("s2.gz")))
        self.assertIsNone(s.cohort)
        self.assertIsNone(s.population)

    def test_metadata_coerced_to_string(self):
        conf = self._base_conf(self._sumstats("s3.gz"))
        conf["cohort"] = 7  # numeric in JSON -> stored as string
        s = Study(conf)
        self.assertEqual(s.cohort, "7")


class TestLeaveOutIntegration(unittest.TestCase):
    """End-to-end runs of meta_analysis.run() with subprocess (bgzip/tabix)
    mocked so the plain output TSV stays on disk for inspection."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def _write_sumstats(self, name):
        """Single shared variant present in every study so they always match."""
        path = os.path.join(self.temp_dir, name)
        with gzip.open(path, "wt") as f:
            f.write("CHR\tPOS\tREF\tALT\tbeta\tpval\n")
            f.write("1\t1000\tA\tG\t0.1\t0.01\n")
        return path

    def _conf_entry(self, name, cohort=None, population=None):
        entry = {
            "name": name, "file": self._write_sumstats(name + ".gz"),
            "n_cases": 100, "n_controls": 200,
            "chr": "CHR", "pos": "POS", "ref": "REF", "alt": "ALT",
            "effect": "beta", "effect_type": "beta", "pval": "pval",
        }
        if cohort is not None:
            entry["cohort"] = cohort
        if population is not None:
            entry["population"] = population
        return entry

    def _run(self, entries, opts):
        conf_path = os.path.join(self.temp_dir, "conf.json")
        with open(conf_path, "w") as f:
            json.dump({"meta": entries}, f)
        out_path = os.path.join(self.temp_dir, "out.tsv")
        argv = ["meta_analysis.py", conf_path, out_path, "n"] + opts
        err = io.StringIO()
        with patch.object(meta_analysis.subprocess, "run"), \
             patch.object(sys, "argv", argv), \
             redirect_stderr(err):
            meta_analysis.run()
        with open(out_path) as f:
            lines = [ln.rstrip("\n") for ln in f]
        header = lines[0].split("\t")
        rows = [ln.split("\t") for ln in lines[1:]]
        return header, rows, err.getvalue()

    def _first_row_dict(self, header, rows):
        self.assertTrue(rows, "expected at least one data row")
        return dict(zip(header, rows[0]))

    def test_no_metadata_skips_grouped_loo(self):
        entries = [self._conf_entry("S1"), self._conf_entry("S2"), self._conf_entry("S3")]
        header, rows, err = self._run(
            entries,
            ["--leave_one_out", "--leave_one_cohort_out", "--leave_one_population_out"],
        )
        # Per-study LOO columns still present
        self.assertIn("leave_S1_n_meta_p", header)
        # No grouped columns
        self.assertFalse([h for h in header if h.startswith("leave_cohort_")])
        self.assertFalse([h for h in header if h.startswith("leave_population_")])
        # Warn-skip messages emitted
        self.assertIn("leave-one-cohort-out", err)
        self.assertIn("leave-one-population-out", err)

    def test_all_metadata_grouped_columns_and_counts(self):
        entries = [
            self._conf_entry("S1", cohort="A", population="EUR"),
            self._conf_entry("S2", cohort="A", population="AFR"),
            self._conf_entry("S3", cohort="B", population="EUR"),
        ]
        header, rows, _ = self._run(
            entries, ["--leave_one_cohort_out", "--leave_one_population_out"]
        )
        for col in ["leave_cohort_A_N", "leave_cohort_B_N",
                    "leave_population_EUR_N", "leave_population_AFR_N"]:
            self.assertIn(col, header)
        d = self._first_row_dict(header, rows)
        # leave cohort A -> drop S1,S2 -> only S3 remains
        self.assertEqual(d["leave_cohort_A_N"], "1")
        # leave cohort B -> drop S3 -> S1,S2 remain
        self.assertEqual(d["leave_cohort_B_N"], "2")
        # leave population EUR -> drop S1,S3 -> only S2 remains
        self.assertEqual(d["leave_population_EUR_N"], "1")
        # leave population AFR -> drop S2 -> S1,S3 remain
        self.assertEqual(d["leave_population_AFR_N"], "2")

    def test_partial_metadata_fieldless_always_retained(self):
        entries = [
            self._conf_entry("S1", cohort="A", population="EUR"),
            self._conf_entry("S2", cohort="A", population="AFR"),
            self._conf_entry("S3", cohort="B", population="EUR"),
            self._conf_entry("S4"),  # no cohort/population -> always retained
        ]
        header, rows, _ = self._run(
            entries, ["--leave_one_cohort_out", "--leave_one_population_out"]
        )
        # No group is created for the field-less study
        self.assertNotIn("leave_cohort_None_N", header)
        self.assertNotIn("leave_population_None_N", header)
        d = self._first_row_dict(header, rows)
        # leave cohort A -> drop S1,S2 -> S3 + field-less S4 remain
        self.assertEqual(d["leave_cohort_A_N"], "2")
        # leave population AFR -> drop S2 -> S1,S3,S4 remain
        self.assertEqual(d["leave_population_AFR_N"], "3")

    def test_backward_compat_study_loo_only(self):
        entries = [self._conf_entry("S1"), self._conf_entry("S2"), self._conf_entry("S3")]
        header, _, _ = self._run(entries, ["--leave_one_out"])
        # Exactly the per-study LOO columns, nothing grouped
        self.assertEqual(
            [h for h in header if h.startswith("leave_")],
            ["leave_S1_N", "leave_S1_n_meta_beta", "leave_S1_n_meta_sebeta",
             "leave_S1_n_meta_p", "leave_S1_n_meta_mlogp",
             "leave_S2_N", "leave_S2_n_meta_beta", "leave_S2_n_meta_sebeta",
             "leave_S2_n_meta_p", "leave_S2_n_meta_mlogp",
             "leave_S3_N", "leave_S3_n_meta_beta", "leave_S3_n_meta_sebeta",
             "leave_S3_n_meta_p", "leave_S3_n_meta_mlogp"],
        )

    def test_no_loo_flags_no_leave_columns(self):
        entries = [self._conf_entry("S1"), self._conf_entry("S2"), self._conf_entry("S3")]
        header, _, _ = self._run(entries, [])
        self.assertFalse([h for h in header if h.startswith("leave_")])


if __name__ == "__main__":
    unittest.main()
