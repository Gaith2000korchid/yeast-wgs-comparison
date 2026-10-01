"""Review the known synthetic SNPs after end-to-end alignment/calling."""
import csv
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
from review_loci import review
truth = json.loads(Path('data/demo/truth.json').read_text())
loci = Path('results/review.synthetic.loci.tsv')
with loci.open('w', newline='') as stream:
    fields = ['chrom','pos_1based','ref','alt','category','illumina_gt','ont_gt']
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for row in truth['snps']:
        writer.writerow({'chrom':row['chrom'],'pos_1based':row['pos'],'ref':row['ref'],'alt':row['alt'],
                         'category':'shared','illumina_gt':'1/1','ont_gt':'1/1'})
review(loci, 'data/demo/reference.fa', 'results/bam/synthetic.illumina.bam',
       'results/bam/synthetic.ont.bam', 'results/review.synthetic')
data=json.loads(Path('results/review.synthetic/read_support.json').read_text())
assert len(data['loci']) == 3
for row in data['loci']:
    for view in row['views'].values():
        assert view['depth'] >= 8, view
        assert view['alt_fraction'] > 0.95, view
assert data['visual_review_status'] == 'not_performed'
print('Read support recovered the three planted SNPs; no visual validation claimed.')
