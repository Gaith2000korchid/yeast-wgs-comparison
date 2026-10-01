"""Plot selected-locus alternate fraction; cells show ALT/retained depth."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def plot(source, output):
    plt.rcParams['svg.hashsalt']='yeast-wgs-review'
    data=json.loads(Path(source).read_text())
    rows=data['loci']
    views=['illumina_bq13','ont_bq13','ont_bq7']
    matrix=np.array([[np.nan if row['views'][v]['alt_fraction'] is None else row['views'][v]['alt_fraction'] for v in views] for row in rows])
    fig,ax=plt.subplots(figsize=(8,8.5))
    cmap=plt.get_cmap('Blues').copy()
    cmap.set_bad('#ededed')
    image=ax.imshow(matrix,aspect='auto',cmap=cmap,vmin=0,vmax=1)
    ax.set_xticks(range(3),['Illumina BQ≥13','ONT BQ≥13','ONT BQ≥7'])
    ax.set_yticks(range(len(rows)),[f"{r['chrom']}:{r['pos_1based']}  {r['ref']}>{r['alt']}" for r in rows])
    for i,row in enumerate(rows):
        for j,v in enumerate(views):
            s=row['views'][v]
            ax.text(j,i,f"{s['alt_count']}/{s['depth']}" if s['depth'] else 'no bases',ha='center',va='center',fontsize=10,color='white' if s['depth'] and s['alt_fraction']>.65 else '#17283b')
    for i in range(1,len(rows)):
        if rows[i]['category'] != rows[i-1]['category']:
            ax.axhline(i-.5,color='white',linewidth=3)
    ax.set_title('P11: original-base support at selected SNPs\nCells: alternate / retained fragments or reads',loc='left',fontsize=12,pad=20)
    colorbar=fig.colorbar(image,ax=ax,fraction=.045,pad=.04)
    colorbar.set_label('Alternate fraction (diagnostic)')
    fig.text(.03,.025,'Groups follow the input locus order (see evidence TSV).\nMAPQ≥20; original qualities, no BAQ; overlapping mates count once.\nSelected convenience sample; visual IGV review pending; no accuracy estimate.',fontsize=9)
    fig.subplots_adjust(left=.27,right=.84,top=.87,bottom=.14)
    fig.savefig(output,dpi=160,metadata={'Date':None})
    plt.close(fig)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',default='results/P11/review/read_support.json')
    parser.add_argument('--output',default='results/P11/review/allele_support.svg')
    a=parser.parse_args()
    plot(a.input,a.output)
