"""Integration assertion: all three toy SNPs must be recovered by both branches."""
import json
from pathlib import Path

truth = json.loads(Path("data/demo/truth.json").read_text())
with open("results/comparison/synthetic.snps.tsv") as stream:
    import csv
    rows = list(csv.DictReader(stream, delimiter="\t"))
expected = {(x["chrom"], x["pos"], x["ref"], x["alt"]) for x in truth["snps"]}
observed = {(x["chrom"], int(x["pos_1based"]), x["ref"], x["alt"]) for x in rows if x["category"] == "shared"}
assert expected <= observed, f"Missing known synthetic SNPs: {expected - observed}"
assert Path("results/report.html").is_file()
assert Path("results/multiqc/multiqc_report.html").is_file()
print("Both branches recovered all 3 known synthetic SNPs; reports exist.")
