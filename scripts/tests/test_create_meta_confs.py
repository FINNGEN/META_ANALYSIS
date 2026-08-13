#!/usr/bin/env python3

import unittest
import sys
import os
import io
import json
import tempfile
import shutil
import subprocess
from types import SimpleNamespace
from contextlib import redirect_stdout

# Add parent directory to path to import the module under test
SCRIPTS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, SCRIPTS_DIR)
import create_meta_confs
from create_meta_confs import get_studies, categorize_columns, generate_json, format_int

SCRIPT_PATH = os.path.join(SCRIPTS_DIR, "create_meta_confs.py")


def col_spec(cohort=None, population=None, se="sebeta"):
    """A per-study column spec as it would appear in --study_meta_json, plus the
    link/n_cases/n_controls mappings that categorize_columns() adds."""
    spec = {
        "chr": "#CHR", "pos": "POS", "ref": "REF", "alt": "ALT",
        "effect": "beta", "effect_type": "beta", "pval": "pval",
        "extra_cols": ["af_alt"],
    }
    if se is not None:
        spec["se"] = se
    if cohort is not None:
        spec["cohort"] = cohort
    if population is not None:
        spec["population"] = population
    return spec


class TestUnitHelpers(unittest.TestCase):

    def test_format_int(self):
        self.assertEqual(format_int("123"), 123)
        self.assertEqual(format_int("NA"), 0)
        self.assertEqual(format_int(""), 0)

    def test_get_studies(self):
        header = ["phenotype", "S1_link", "S1_n_cases", "S2_link", "S2_n_controls"]
        self.assertEqual(get_studies(header, "_link"), ["S1", "S2"])

    def test_get_studies_no_match(self):
        header = ["phenotype", "S1_n_cases", "S2_n_controls"]
        self.assertEqual(get_studies(header, "_link"), [])

    def test_categorize_columns_keeps_metadata(self):
        # cohort/population provided in the spec json must survive categorization
        meta = {"S1": col_spec(cohort="FG", population="EUR"), "S2": col_spec()}
        tmp = tempfile.mkdtemp()
        try:
            spec_path = os.path.join(tmp, "spec.json")
            with open(spec_path, "w") as f:
                json.dump(meta, f)
            header = ["phenotype",
                      "S1_link", "S1_n_cases", "S1_n_controls",
                      "S2_link", "S2_n_cases", "S2_n_controls"]
            args = SimpleNamespace(
                study_meta_json=spec_path, studies=["S1", "S2"],
                phenotype_col="phenotype", link_col_suffix="_link",
                n_cases_col_suffix="_n_cases", n_controls_col_suffix="_n_controls",
            )
            columns = categorize_columns(header, args)
            self.assertEqual(columns["phenotype"], "phenotype")
            self.assertEqual(columns["studies"]["S1"]["link"], "S1_link")
            self.assertEqual(columns["studies"]["S1"]["n_cases"], "S1_n_cases")
            self.assertEqual(columns["studies"]["S1"]["n_controls"], "S1_n_controls")
            self.assertEqual(columns["studies"]["S1"]["cohort"], "FG")
            self.assertEqual(columns["studies"]["S1"]["population"], "EUR")
            self.assertNotIn("cohort", columns["studies"]["S2"])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


class TestGenerateJson(unittest.TestCase):
    """generate_json writes <phenotype>.json into JSONS_DIR if it exists, else
    the cwd. Run inside a tempdir (no jsons/ dir) and read the file back."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.cwd = os.getcwd()
        os.chdir(self.tmp)
        self.header = ["phenotype",
                       "S1_link", "S1_n_cases", "S1_n_controls",
                       "S2_link", "S2_n_cases", "S2_n_controls"]
        self.h_idx = {h: i for i, h in enumerate(self.header)}
        self.data = ["pheno1",
                     "gs://bucket/S1.gz", "100", "200",
                     "gs://bucket/S2.gz", "300", "400"]

    def tearDown(self):
        os.chdir(self.cwd)
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _columns(self, s1_spec, s2_spec):
        cols = {"phenotype": "phenotype", "studies": {"S1": dict(s1_spec), "S2": dict(s2_spec)}}
        for s in ("S1", "S2"):
            cols["studies"][s]["link"] = s + "_link"
            cols["studies"][s]["n_cases"] = s + "_n_cases"
            cols["studies"][s]["n_controls"] = s + "_n_controls"
        return cols

    def _run(self, columns, continuous=False):
        generate_json(self.data, self.h_idx, columns, ["S1", "S2"], continuous)
        with open("pheno1.json") as f:
            return json.load(f)["meta"]

    def test_core_fields_and_gcs_rewrite(self):
        meta = self._run(self._columns(col_spec(), col_spec()))
        s1 = meta[0]
        self.assertEqual(s1["name"], "S1")
        self.assertEqual(s1["file"], "/cromwell_root/bucket/S1.gz")
        self.assertEqual(s1["n_cases"], 100)
        self.assertEqual(s1["n_controls"], 200)
        self.assertEqual(s1["chr"], "#CHR")
        self.assertEqual(s1["pval"], "pval")
        self.assertEqual(s1["se"], "sebeta")
        self.assertEqual(s1["extra_cols"], ["af_alt"])

    def test_cohort_population_present(self):
        meta = self._run(self._columns(col_spec(cohort="FG", population="EUR"), col_spec()))
        self.assertEqual(meta[0]["cohort"], "FG")
        self.assertEqual(meta[0]["population"], "EUR")

    def test_cohort_population_absent_omitted(self):
        meta = self._run(self._columns(col_spec(), col_spec()))
        self.assertNotIn("cohort", meta[0])
        self.assertNotIn("population", meta[0])

    def test_partial_metadata(self):
        # Only S1 tagged; S2 left untagged
        meta = self._run(self._columns(col_spec(cohort="FG", population="EUR"), col_spec()))
        self.assertIn("cohort", meta[0])
        self.assertIn("population", meta[0])
        self.assertNotIn("cohort", meta[1])
        self.assertNotIn("population", meta[1])

    def test_se_optional(self):
        meta = self._run(self._columns(col_spec(se=None), col_spec()))
        self.assertNotIn("se", meta[0])

    def test_continuous_zero_controls(self):
        meta = self._run(self._columns(col_spec(), col_spec()), continuous=True)
        self.assertEqual(meta[0]["n_controls"], 0)
        self.assertEqual(meta[1]["n_controls"], 0)


class TestIntegration(unittest.TestCase):
    """End-to-end runs of the script via subprocess in a tempdir."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _write_mapping(self, rows):
        path = os.path.join(self.tmp, "mapping.tsv")
        header = ["phenotype",
                  "S1_link", "S1_n_cases", "S1_n_controls",
                  "S2_link", "S2_n_cases", "S2_n_controls",
                  "S3_link", "S3_n_cases", "S3_n_controls"]
        with open(path, "w") as f:
            f.write("\t".join(header) + "\n")
            for r in rows:
                f.write("\t".join(r) + "\n")
        return path

    def _write_spec(self):
        path = os.path.join(self.tmp, "spec.json")
        spec = {
            "S1": col_spec(cohort="FG", population="EUR"),
            "S2": col_spec(cohort="UKBB", population="AFR"),
            "S3": col_spec(),  # untagged
        }
        with open(path, "w") as f:
            json.dump(spec, f)
        return path

    def _run(self, mapping, extra_args):
        sums = os.path.join(self.tmp, "sums.txt")
        jsons = os.path.join(self.tmp, "jsons.txt")
        cmd = [sys.executable, SCRIPT_PATH, mapping,
               "--sumstat_filelist_name", sums,
               "--json_filelist_name", jsons] + extra_args
        return subprocess.run(cmd, cwd=self.tmp, capture_output=True, text=True)

    def _row(self, pheno, s1="gs://b/S1.gz", s2="gs://b/S2.gz", s3="gs://b/S3.gz"):
        return [pheno, s1, "100", "200", s2, "300", "400", s3, "500", "600"]

    def test_passthrough_end_to_end(self):
        mapping = self._write_mapping([self._row("phenoA")])
        spec = self._write_spec()
        res = self._run(mapping, ["--study_meta_json", spec])
        self.assertEqual(res.returncode, 0, res.stderr)
        out = os.path.join(self.tmp, "jsons", "phenoA.json")
        self.assertTrue(os.path.exists(out), res.stdout + res.stderr)
        with open(out) as f:
            meta = {s["name"]: s for s in json.load(f)["meta"]}
        self.assertEqual(meta["S1"]["cohort"], "FG")
        self.assertEqual(meta["S1"]["population"], "EUR")
        self.assertEqual(meta["S2"]["cohort"], "UKBB")
        self.assertNotIn("cohort", meta["S3"])
        self.assertEqual(meta["S1"]["file"], "/cromwell_root/b/S1.gz")

    def test_min_studies_skips_phenotype(self):
        # phenoB has only S1 present -> below default min_studies (2) -> skipped
        mapping = self._write_mapping([self._row("phenoB", s2="NA", s3="NA")])
        spec = self._write_spec()
        res = self._run(mapping, ["--study_meta_json", spec])
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "jsons", "phenoB.json")))

    def test_required_studies_skips_phenotype(self):
        # phenoC lacks S1; require S1 -> skipped
        mapping = self._write_mapping([self._row("phenoC", s1="NA")])
        spec = self._write_spec()
        res = self._run(mapping, ["--study_meta_json", spec, "--required_studies", "S1"])
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "jsons", "phenoC.json")))

    def test_complete_requires_all_studies(self):
        # --complete: phenoD missing S3 -> skipped; phenoE complete -> kept
        mapping = self._write_mapping([
            self._row("phenoD", s3="NA"),
            self._row("phenoE"),
        ])
        spec = self._write_spec()
        res = self._run(mapping, ["--study_meta_json", spec, "--complete"])
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "jsons", "phenoD.json")))
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "jsons", "phenoE.json")))

    def test_continuous_flag(self):
        mapping = self._write_mapping([self._row("phenoF")])
        spec = self._write_spec()
        res = self._run(mapping, ["--study_meta_json", spec, "--continuous"])
        self.assertEqual(res.returncode, 0, res.stderr)
        with open(os.path.join(self.tmp, "jsons", "phenoF.json")) as f:
            for s in json.load(f)["meta"]:
                self.assertEqual(s["n_controls"], 0)

    def test_old_flag_rejected(self):
        mapping = self._write_mapping([self._row("phenoA")])
        spec = self._write_spec()
        res = self._run(mapping, ["--study_cols_json", spec])
        self.assertNotEqual(res.returncode, 0)
        self.assertIn("unrecognized arguments", res.stderr)


if __name__ == "__main__":
    unittest.main()
