"""Export compact measured P11 QC/overlap evidence for the GitHub repository.

Requires matplotlib only for plotting; the workflow itself does not use it.
Run after the entire real workflow has completed successfully.
"""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import zipfile


def load(path):
    return json.loads(Path(path).read_text())


def nuclear(coverage):
    chromosomes = {k:v for k,v in coverage['by_chromosome'].items() if k != 'Mito'}
    bases = sum(v['bases'] for v in chromosomes.values())
    return {'reference_bases': bases,
            'mean_filtered_depth': sum(v['mean_depth']*v['bases'] for v in chromosomes.values())/bases,
            'covered_fraction': sum(v['covered_bases'] for v in chromosomes.values())/bases,
            'depth_eligible_fraction': sum(v['depth_eligible_bases'] for v in chromosomes.values())/bases}


def alignment(path):
    text=Path(path).read_text()
    values={}
    for key,label in [('primary_reads','primary'),('primary_mapped','primary mapped'),('primary_duplicates','primary duplicates')]:
        match=re.search(r'^(\d+) \+ (\d+) '+label+r'(?: \(|$)',text,re.M)
        if not match: raise ValueError(f'Missing flagstat metric: {label}')
        values[key]=sum(map(int,match.groups()))
    values['primary_mapped_fraction']=values['primary_mapped']/values['primary_reads']
    values['primary_duplicate_fraction']=values['primary_duplicates']/values['primary_reads']
    return values


def allele_balance(path):
    count=het=hom_alt=0
    bins=[0]*20
    with gzip.open(path,'rt') as f:
        for line in f:
            if line.startswith('#'):continue
            chrom,pos,_,ref,alt,qual,filt,info,fmt,sample=line.rstrip().split('\t')
            if chrom=='Mito' or len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT':continue
            values=dict(zip(fmt.split(':'),sample.split(':')))
            gt=values.get('GT','.').replace('|','/').split('/')
            if '.' in gt or '1' not in gt:continue
            count+=1
            if sorted(gt)==['0','1']:
                het+=1
                ad=values.get('AD','.').split(',')
                if len(ad)==2 and all(a.isdigit() for a in ad) and sum(map(int,ad))>0:
                    balance=int(ad[1])/sum(map(int,ad))
                    bins[min(int(balance*20),19)]+=1
            elif gt==['1','1']:hom_alt+=1
    return {'nuclear_alternate_snp_records':count,'heterozygous_snp_records':het,'homozygous_alt_snp_records':hom_alt,
            'heterozygous_alt_balance_bins':bins,'bin_width':0.05}


def export(results,prepared,out):
    root=Path(results);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if not (root/'report.html').is_file():raise ValueError('Complete the workflow before exporting')
    prep=load(Path(prepared)/'preparation.json')
    cover={p:load(root/f'qc/P11/{p}.coverage.json') for p in ['illumina','ont']}
    comparison=load(root/'comparison/P11.comparison.json')
    fastp=load(root/'qc/P11/fastp.json')['summary']
    metrics={'sample':'P11 / CIC / Ponton11','sample_accession':'SAMEA3895683','execution_date':'2026-10-01',
             'ploidy':2,'illumina_preparation':prep['illumina'],
             'ont_selected_2d_reads':sum(v['selected_2d_reads'] for v in prep['ont_runs']),
             'ont_selected_2d_bases':sum(v['selected_2d_bases'] for v in prep['ont_runs']),
             'fastp_summary':fastp,'nuclear_coverage':{p:nuclear(cover[p]) for p in cover},
             'alignment':{p:alignment(root/f'qc/P11/{p}.flagstat.txt') for p in cover},
             'nuclear_snp_genotypes':{p:allele_balance(root/f'variants/P11.{p}.pass.vcf.gz') for p in cover},
             'comparison':comparison}
    metrics['shared_nuclear_depth_fraction']=comparison['shared_depth_eligible_bases']/metrics['nuclear_coverage']['illumina']['reference_bases']
    (out/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    shutil.copyfile(Path(prepared)/'preparation.json',out/'preparation.json')
    wanted=[root/'comparison/P11.comparison.json',root/'provenance/tool_versions.txt']
    for platform in cover:
        wanted.extend(root/f'qc/P11/{platform}.{suffix}' for suffix in ['coverage.json','flagstat.txt','samtools.stats.txt'])
        wanted.append(root/f'variants/P11.{platform}.bcftools.stats.txt')
    wanted.extend((root/'qc/P11/nanoplot').glob('*Stats*.txt'))
    for src in wanted:shutil.copyfile(src,out/src.name)
    # Preserve all raw/clean FastQC module outcomes in a compact reviewable table.
    with (out/'fastqc_checks.tsv').open('w') as dest:
        dest.write('stage\tfile\tstatus\tmodule\n')
        for path in sorted((root/'qc/P11/fastqc').glob('*/*.zip')):
            with zipfile.ZipFile(path) as z:
                name=next(n for n in z.namelist() if n.endswith('/summary.txt'))
                for line in z.read(name).decode().splitlines():
                    status,module,filename=line.split('\t')
                    dest.write(f'{path.parent.name}\t{filename}\t{status}\t{module}\n')
    chromosomes=list(cover['illumina']['by_chromosome'])
    with (out/'coverage_by_chromosome.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['chromosome','platform','bases','mean_filtered_depth','depth_eligible_fraction'])
        for chrom in chromosomes:
            for platform in cover:
                v=cover[platform]['by_chromosome'][chrom]
                w.writerow([chrom,platform,v['bases'],v['mean_depth'],v['depth_eligible_bases']/v['bases']])
    # Illustrative loci, selected deterministically; no claim that they are errors.
    examples=[];counts={}
    with (root/'comparison/P11.snps.tsv').open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            cat=row['category']
            if counts.get(cat,0)<5:
                examples.append(row);counts[cat]=counts.get(cat,0)+1
    with (out/'loci_for_review.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=['chrom','pos_1based','ref','alt','category','illumina_gt','ont_gt'],delimiter='\t')
        w.writeheader();w.writerows(examples)
    # Plot all nuclear contigs in reference order, explicitly excluding Mito.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'svg.hashsalt':'yeast-wgs-P11'})
    chroms=[c for c in chromosomes if c!='Mito'];x=list(range(len(chroms)))
    fig,axes=plt.subplots(2,1,figsize=(10,6),sharex=True,layout='constrained')
    for platform,color,offset in [('illumina','#245f99',-.18),('ont','#d76b28',.18)]:
        values=cover[platform]['by_chromosome']
        axes[0].bar([v+offset for v in x],[values[c]['mean_depth'] for c in chroms],width=.35,label=platform,color=color)
        axes[1].bar([v+offset for v in x],[100*values[c]['depth_eligible_bases']/values[c]['bases'] for c in chroms],width=.35,color=color)
    axes[0].set_ylabel('Mean filtered depth (×)');axes[0].legend(frameon=False)
    axes[1].set_ylabel('Depth eligible (%)');axes[1].set_ylim(0,105)
    axes[1].set_xticks(x,chroms);axes[1].set_xlabel('Nuclear chromosome')
    for ax in axes:ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
    fig.suptitle('CIC / Ponton11: Illumina 25% pairs vs historical ONT 2D\nMAPQ ≥20 · base quality ≥13 · eligible depth 8–200',fontsize=12)
    fig.savefig(out/'coverage.svg',metadata={'Date':None})
    fig.savefig(root/'coverage_preview.png',dpi=150)
    plt.close(fig)
    with (out/'SHA256SUMS').open('w') as f:
        for path in sorted(out.iterdir()):
            if path.name!='SHA256SUMS' and path.is_file():
                f.write(hashlib.sha256(path.read_bytes()).hexdigest()+'  '+path.name+'\n')
    print(json.dumps(metrics,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',default='results/P11')
    p.add_argument('--prepared',default='data/real/prepared')
    p.add_argument('--out',default='docs/results/P11')
    a=p.parse_args();export(a.results,a.prepared,a.out)
