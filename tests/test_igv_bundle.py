"""Regional export must preserve SNP evidence and avoid duplicate records."""
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
import zipfile
import pysam
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow/scripts'))
from export_igv_bundle import export, windows


class BundleTests(unittest.TestCase):
    def test_portable_export_preserves_evidence_and_original_records(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            ref=root/'ref.fa'
            ref.write_text('>I\n'+'A'*100+'\n')
            pysam.faidx(str(ref))
            loci=root/'loci.tsv'
            loci.write_text('chrom\tpos_1based\tref\talt\tcategory\tillumina_gt\tont_gt\nI\t10\tA\tG\tshared\t1/1\t1/1\nI\t20\tA\tG\tshared\t1/1\t1/1\n')
            bam=root/'source.bam'
            with pysam.AlignmentFile(str(bam),'wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'I','LN':100}]}) as stream:
                for name,start,length in [('spans_both',0,30),('outside',80,10)]:
                    read=pysam.AlignedSegment()
                    read.query_name=name
                    seq=list('A'*length)
                    if start==0: seq[9]=seq[19]='G'
                    read.query_sequence=''.join(seq)
                    read.query_qualities=[20]*length
                    read.reference_id=0;read.reference_start=start;read.mapping_quality=60;read.cigarstring=f'{length}M'
                    stream.write(read)
            pysam.index(str(bam))
            archive=export(loci,ref,bam,bam,root/'bundle',flank=0,source_commit='fixture')
            manifest=json.loads((root/'bundle/manifest.json').read_text())
            self.assertEqual(manifest['alignment_record_counts'],{'illumina':1,'ont':1})
            self.assertTrue(manifest['selected_locus_evidence_matches_source'])
            self.assertEqual(manifest['visual_review_status'],'not_performed')
            with pysam.AlignmentFile(str(root/'bundle/P11.ont.regional.bam'),'rb') as stream:
                record=next(stream.fetch())
                self.assertEqual(record.cigarstring,'30M')
                self.assertEqual(record.query_length,30)
            relocated=root/'different folder'
            with zipfile.ZipFile(archive) as zipped: zipped.extractall(relocated)
            for line in (relocated/'bundle/SHA256SUMS').read_text().splitlines():
                expected,path=line.split('  ',1)
                self.assertEqual(hashlib.sha256((relocated/'bundle'/path).read_bytes()).hexdigest(),expected)
            with self.assertRaisesRegex(ValueError,'already exist'):
                export(loci,ref,bam,bam,root/'bundle')
            with self.assertRaisesRegex(ValueError,'nonnegative'):
                windows([],ref,-1)
            self.assertEqual(windows([{'chrom':'I','pos_1based':1,'ref':'A'},{'chrom':'I','pos_1based':5,'ref':'A'}],ref,5),[('I',0,10)])

if __name__=='__main__':unittest.main()
