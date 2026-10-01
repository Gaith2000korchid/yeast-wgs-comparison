"""Regression smoke test: NanoPlot must plot more than its 10,000-point limit.

A 400-read workflow toy fixture did not expose the Plotly >=6 length mismatch
in historical NanoPlot 1.43. This fixture is software test data only.
"""
import gzip
from pathlib import Path
import subprocess

out=Path('results/nanoplot_regression')
out.mkdir(parents=True,exist_ok=True)
reads=out/'fixture.fastq.gz'
with gzip.open(reads,'wt') as f:
    for i in range(12001):
        length=100+i%600
        seq=('ACGT'*175)[:length]
        qual=chr(33+8+i%20)*length
        f.write(f'@qc_fixture_{i}\n{seq}\n+\n{qual}\n')
with (out/'execution.log').open('w') as log:
    subprocess.run(['NanoPlot','--fastq',str(reads),'--threads','2','--N50','--outdir',str(out/'qc')],stdout=log,stderr=subprocess.STDOUT,check=True)
stats=(out/'qc/NanoStats.txt').read_text()
assert '12,001' in stats or '12001' in stats,stats
assert (out/'qc/NanoPlot-report.html').is_file()
print('NanoPlot >10,000-read regression passed (12,001 synthetic reads)')
