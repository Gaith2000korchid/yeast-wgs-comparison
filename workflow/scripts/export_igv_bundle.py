"""Export indexed regional BAMs and a complete reference for portable IGV review.

Original records overlapping merged locus windows are retained once with the
multi-region iterator. Evidence at selected sites must match the source BAMs.
The regional files cannot be used for genome-wide QC or variant calling.
"""
import argparse
import json
from pathlib import Path
import shutil
import tempfile
import xml.etree.ElementTree as ET
import zipfile
import pysam
from review_loci import digest, igv_session, read_loci, review


def windows(rows, reference, flank):
    if flank < 0:
        raise ValueError('Flank must be nonnegative')
    intervals = {}
    with pysam.FastaFile(str(reference)) as fasta:
        for row in rows:
            chrom, pos = row['chrom'], row['pos_1based'] - 1
            if chrom not in fasta.references or pos >= fasta.get_reference_length(chrom):
                raise ValueError('Locus outside reference')
            if fasta.fetch(chrom,pos,pos+1).upper() != row['ref']:
                raise ValueError('Reference allele mismatch')
            intervals.setdefault(chrom,[]).append((max(0,pos-flank),min(fasta.get_reference_length(chrom),pos+flank+1)))
        merged = []
        for chrom in fasta.references:
            blocks = []
            for start,end in sorted(intervals.get(chrom,[])):
                if blocks and start <= blocks[-1][1]:
                    blocks[-1] = (blocks[-1][0],max(blocks[-1][1],end))
                else:
                    blocks.append((start,end))
            merged.extend((chrom,start,end) for start,end in blocks)
    return merged


def export(loci, reference, illumina, ont, output, flank=2000, source_commit='unknown'):
    loci,reference,illumina,ont=map(Path,(loci,reference,illumina,ont))
    out=Path(output)
    archive=Path(str(out)+'.zip')
    if out.exists() or archive.exists():
        raise ValueError('Output directory must not already exist')
    rows=read_loci(loci)
    regions=windows(rows,reference,flank)
    with tempfile.TemporaryDirectory() as temp:
        before=Path(temp)/'source_evidence'
        review(loci,reference,illumina,ont,before)
        source=json.loads((before/'read_support.json').read_text())
    out.mkdir(parents=True)
    shutil.copyfile(reference,out/'reference.fa')
    pysam.faidx(str(out/'reference.fa'))
    shutil.copyfile(loci,out/'loci.tsv')
    bed=out/'extraction_windows.bed'
    bed.write_text(''.join(f'{c}\t{s}\t{e}\n' for c,s,e in regions))
    counts={}
    for platform,bam in [('illumina',illumina),('ont',ont)]:
        target=out/f'P11.{platform}.regional.bam'
        pysam.view('-b','-M','-L',str(bed),'-o',str(target),str(bam),catch_stdout=False)
        pysam.index(str(target))
        pysam.quickcheck(str(target))
        with pysam.AlignmentFile(str(target),'rb') as stream:
            counts[platform]=sum(1 for _ in stream.fetch(until_eof=True))
    review(out/'loci.tsv',out/'reference.fa',out/'P11.illumina.regional.bam',out/'P11.ont.regional.bam',out/'evidence')
    after=json.loads((out/'evidence/read_support.json').read_text())
    if source['loci'] != after['loci']:
        raise ValueError('Regional BAM evidence differs from source BAM evidence')
    focus=next((r for r in rows if r['chrom']=='III' and r['pos_1based']==105658),rows[0])
    with pysam.FastaFile(str(out/'reference.fa')) as fasta:
        length=fasta.get_reference_length(focus['chrom'])
    locus=f"{focus['chrom']}:{max(1,focus['pos_1based']-100)}-{min(length,focus['pos_1based']+100)}"
    igv_session(out/'session.xml',out/'reference.fa',{'Illumina regional':out/'P11.illumina.regional.bam','ONT historical 2D regional':out/'P11.ont.regional.bam'},locus)
    xml=ET.parse(out/'session.xml').getroot()
    for path in [xml.attrib['genome']]+[r.attrib[k] for r in xml.findall('Resources/Resource') for k in ('path','index')]:
        if not (out/path).resolve().is_file():
            raise ValueError('IGV resource path does not resolve')
    manifest={'source_commit':source_commit,'pysam_version':pysam.__version__,'script_sha256':digest(__file__),
              'flank_bp':flank,'loci':len(rows),'merged_windows':regions,'alignment_record_counts':counts,
              'source_fingerprints':source['sources'],'selected_locus_evidence_matches_source':True,
              'visual_review_status':'not_performed','scope':'Only records overlapping the extraction windows; full original alignments retained once. Coverage outside windows is incomplete. No genome-wide inference or re-calling from this bundle.'}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (out/'LIRE_MOI.txt').write_text(
        'P11 - Inspection IGV, pas de validation visuelle encore effectuee.\n\n'
        '1. Extraire tout le ZIP dans un dossier (ne pas ouvrir depuis le ZIP).\n'
        '2. Dans IGV Desktop : File > Open Session > session.xml.\n'
        '3. Examiner III:105558-105758 puis IV:308149-308349.\n'
        '4. Garder les BAM, BAI, reference.fa et FAI dans ce dossier.\n'
        '5. Renseigner evidence/visual_review_template.tsv et sauvegarder les captures.\n\n'
        f'Les BAM sont regionaux, avec des fenetres de +/- {flank} bases.\n'
        'Les alignements originaux sont conserves entiers : hors fenetres, la couverture est incomplete.\n'
        'Ne pas utiliser pour une QC globale, un appel de variants ou une estimation de precision.\n'
        'Les comptages diagnostiques ne sont pas les comptages affiches par IGV.\n'
        'Guide complet : https://github.com/Gaith2000korchid/yeast-wgs-comparison/blob/main/docs/IGV_REVIEW_FR.md\n')
    (out/'SHA256SUMS').write_text(''.join(f'{digest(p)}  {p.relative_to(out).as_posix()}\n' for p in sorted(out.rglob('*')) if p.is_file() and p.name!='SHA256SUMS'))
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as zipped:
        for p in sorted(out.rglob('*')):
            if p.is_file():
                zipped.write(p,(Path(out.name)/p.relative_to(out)).as_posix())
    print(f'Portable IGV bundle: {archive}; source-locus evidence matches; visual review pending.')
    return archive


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--loci',default='docs/results/P11/bq7.loci_for_review.tsv')
    parser.add_argument('--reference',default='data/real/reference.fa')
    parser.add_argument('--illumina',default='results/P11/bam/P11.illumina.bam')
    parser.add_argument('--ont',default='results/P11/bam/P11.ont.bam')
    parser.add_argument('--out',default='results/P11/igv_bundle')
    parser.add_argument('--flank',type=int,default=2000)
    parser.add_argument('--source-commit',default='unknown')
    a=parser.parse_args()
    export(a.loci,a.reference,a.illumina,a.ont,a.out,a.flank,a.source_commit)
