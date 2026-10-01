"""On toy Q40 reads, BQ13 and BQ7 must yield identical comparison metrics."""
import json
from pathlib import Path
primary=json.loads(Path('results/comparison/synthetic.comparison.json').read_text())
sensitivity=json.loads(Path('results/sensitivity/ont_bq7/synthetic.comparison.json').read_text())
# Allow a cached primary report from before explicit excluded-contig metadata.
primary.setdefault("excluded_contigs",[])
assert primary==sensitivity,(primary,sensitivity)
print('ONT base-quality sensitivity integration passed on synthetic Q40 fixture')
