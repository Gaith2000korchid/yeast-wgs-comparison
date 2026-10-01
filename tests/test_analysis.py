import gzip
import hashlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts"))
from compare_snps import compare, intersect_masks, read_bed
from coverage import summarize
from generate_demo import generate
from validate_inputs import load_samples


def vcf(path, records):
    header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    path.write_text(header + "".join(f"chr1\t{p}\t.\t{ref}\t{alt}\t60\t{filt}\t.\tGT\t{gt}\n"
                                    for p, ref, alt, gt, filt in records))


class AnalysisTests(unittest.TestCase):
    def test_depth_mask_coordinates_zeros_and_high_depth(self):
        bed = io.StringIO()
        result = summarize(io.StringIO("chr1\t1\t0\nchr1\t2\t8\nchr1\t3\t9\nchr1\t4\t201\nchr2\t1\t10\n"), bed, 8, 200)
        self.assertEqual(bed.getvalue(), "chr1\t1\t3\nchr2\t0\t1\n")
        self.assertEqual(result["depth_eligible_bases"], 3)
        self.assertEqual(result["reference_bases"], 5)
        self.assertEqual(result["covered_fraction"], 0.8)
        self.assertEqual(result["mean_depth"], 228 / 5)

    def test_depth_rejects_missing_positions(self):
        with self.assertRaises(ValueError):
            summarize(io.StringIO("chr1\t2\t8\n"), io.StringIO(), 8, 200)

    def test_mask_intersection(self):
        self.assertEqual(intersect_masks({"chr1": [(0, 5), (10, 15)]}, {"chr1": [(3, 12)]}),
                         {"chr1": [(3, 5), (10, 12)]})

    def test_snps_require_alternate_genotype_and_shared_depth(self):
        with tempfile.TemporaryDirectory() as tmp:
            a, b = Path(tmp) / "a.vcf", Path(tmp) / "b.vcf"
            vcf(a, [(1, "A", "C", "1/1", "PASS"), (2, "G", "T", "0/0", "PASS"),
                    (3, "A", "G", "0/1", "PASS"), (4, "C", "T", "1/1", "LowQual"),
                    (5, "A", "AT", "1/1", "PASS"), (6, "C", "G", "./1", "PASS"),
                    (10, "G", "T", "1/1", "PASS")])
            vcf(b, [(1, "A", "C", "1/1", "PASS"), (3, "A", "G", "1/1", "PASS"),
                    (4, "C", "G", "1/1", "PASS")])
            result, _, _ = compare(a, b, {"chr1": [(0, 4)]})
            self.assertEqual(result["illumina_snps_in_shared_mask"], 2)
            self.assertEqual(result["illumina_alternate_snps_all_regions"], 3)
            self.assertEqual(result["shared_snps"], 2)
            self.assertEqual(result["ont_only_snps"], 1)
            self.assertAlmostEqual(result["snp_jaccard"], 2 / 3)
            self.assertEqual(result["genotype_agreement_on_shared_snps"], 0.5)

    def test_excluded_contig_is_absent_from_domain_and_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "a.vcf"
            vcf(path, [(1, "A", "C", "0/1", "PASS")])
            result, left, right = compare(path, path, {"chr1": [(0, 10)]}, ["chr1"])
            self.assertEqual(result["shared_depth_eligible_bases"], 0)
            self.assertEqual(result["illumina_alternate_snps_all_regions"], 0)
            self.assertEqual(left, {})
            self.assertEqual(right, {})

    def test_empty_overlap_is_not_perfect_agreement(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "a.vcf"
            vcf(path, [])
            result, _, _ = compare(path, path, {})
            self.assertIsNone(result["snp_jaccard"])
            self.assertIsNone(result["genotype_agreement_on_shared_snps"])

    def test_bed_boundary_at_last_included_base(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "a.vcf"
            vcf(path, [(3, "A", "C", "1", "PASS"), (4, "A", "C", "1", "PASS")])
            result, _, _ = compare(path, path, {"chr1": [(2, 3)]})
            self.assertEqual(result["shared_snps"], 1)

    def test_bed_rejects_overlaps(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "a.bed"
            path.write_text("chr1\t0\t10\nchr1\t5\t12\n")
            with self.assertRaises(ValueError):
                read_bed(path)

    def test_sample_sheet_rejects_unmatched_ploidy_and_missing_platform(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "samples.tsv"
            header = "sample\tplatform\tread1\tread2\tploidy\n"
            path.write_text(header + "A\tillumina\tR1\tR2\t1\n")
            with self.assertRaises(ValueError):
                load_samples(path)
            path.write_text(header + "A\tillumina\tR1\tR2\t1\nA\tont\tONT\t\t2\n")
            with self.assertRaises(ValueError):
                load_samples(path)

    def test_demo_is_deterministic_and_has_three_known_snps(self):
        with tempfile.TemporaryDirectory() as tmp:
            left, right = Path(tmp) / "a", Path(tmp) / "b"
            generate(left)
            generate(right)
            for filename in ("reference.fa", "illumina_R1.fastq.gz", "illumina_R2.fastq.gz", "ont.fastq.gz", "truth.json"):
                self.assertEqual(hashlib.sha256((left / filename).read_bytes()).digest(),
                                 hashlib.sha256((right / filename).read_bytes()).digest())
            truth = json.loads((left / "truth.json").read_text())
            self.assertEqual(len(truth["snps"]), 3)
            with gzip.open(left / "illumina_R1.fastq.gz", "rt") as stream:
                lines = stream.readlines()
            self.assertEqual(len(lines), 3000 * 4)
            self.assertEqual(len(lines[1].strip()), len(lines[3].strip()))


if __name__ == "__main__":
    unittest.main()
