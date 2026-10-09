import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/scripts"))
from prepare_local_workflow import prepare


class LocalWorkflowTests(unittest.TestCase):
    def test_only_benchmark_directives_are_removed_and_hashed(self):
        with tempfile.TemporaryDirectory() as folder:
            source, output = Path(folder) / "Snakefile", Path(folder) / "local.smk"
            prefix = b'rule map:\n    input: "reads.fq"\n'
            suffix = b'    shell: "mapper {input:q}"\n'
            source.write_bytes(prefix + b'    benchmark: f"{BENCH}/map.tsv"\n' + suffix)
            record = prepare(source, output)
            self.assertEqual(output.read_bytes(), prefix + suffix)
            self.assertEqual(record["source_sha256"], hashlib.sha256(source.read_bytes()).hexdigest())
            self.assertEqual(record["output_sha256"], hashlib.sha256(output.read_bytes()).hexdigest())
            self.assertEqual(record["removed_benchmark_lines"], [3])
            self.assertEqual(json.loads(output.with_suffix(".smk.provenance.json").read_text()), record)

    def test_multiline_benchmark_is_rejected_without_creating_output(self):
        with tempfile.TemporaryDirectory() as folder:
            source, output = Path(folder) / "Snakefile", Path(folder) / "local.smk"
            source.write_text('rule map:\n    benchmark:\n        "map.tsv"\n')
            with self.assertRaisesRegex(ValueError, "single-line"):
                prepare(source, output)
            self.assertFalse(output.exists())

    def test_source_and_relative_path_context_are_preserved(self):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / "Snakefile"
            source.write_text("rule all:\n    input: []\n")
            with self.assertRaisesRegex(ValueError, "different"):
                prepare(source, source)
            with self.assertRaisesRegex(ValueError, "beside"):
                prepare(source, Path(folder) / "subdir/local.smk")
