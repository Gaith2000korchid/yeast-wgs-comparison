"""Run the portable exporter on both synthetic workflow BAMs."""
from pathlib import Path
import json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow/scripts'))
from export_igv_bundle import export
archive=export('results/review.synthetic.loci.tsv','data/demo/reference.fa',
               'results/bam/synthetic.illumina.bam','results/bam/synthetic.ont.bam',
               'results/igv.synthetic',source_commit='synthetic integration fixture')
manifest=json.loads(Path('results/igv.synthetic/manifest.json').read_text())
assert archive.is_file()
assert manifest['loci']==3
assert manifest['selected_locus_evidence_matches_source']
assert manifest['visual_review_status']=='not_performed'
print('Portable ZIP reproduces full-BAM evidence at all three synthetic SNPs.')
