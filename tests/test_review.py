"""Diagnostic evidence tests using real indexed synthetic BAMs."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
import pysam
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
from review_loci import collect, observation, support, review, read_loci, igv_session


def aligned(name, base='G', quality=20, flag=0, mq=60, cigar='20M'):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = 'A' * 9 + base + 'A' * 10
    read.query_qualities = [quality] * 20
    read.reference_id = 0
    read.reference_start = 0
    read.flag = flag
    read.mapping_quality = mq
    read.cigarstring = cigar
    return read


class ReviewTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.ref = self.root / 'reference.fa'
        self.ref.write_text('>I\n' + 'A' * 100 + '\n')
        pysam.faidx(str(self.ref))
        self.loci = self.root / 'loci.tsv'
        self.loci.write_text('chrom\tpos_1based\tref\talt\tcategory\tillumina_gt\tont_gt\nI\t10\tA\tG\tshared\t0/1\t1/1\n')

    def bam(self, name, reads):
        path = self.root / (name + '.bam')
        with pysam.AlignmentFile(str(path), 'wb', header={'HD': {'VN':'1.6','SO':'coordinate'}, 'SQ':[{'SN':'I','LN':100}]}) as stream:
            for read in reads:
                stream.write(read)
        pysam.index(str(path))
        return path

    def test_fragment_overlap_conflict_quality_and_flags(self):
        reads = [aligned('agree', flag=65), aligned('agree', flag=129),
                 aligned('conflict', flag=65), aligned('conflict', base='A', flag=129),
                 aligned('lowbq', quality=7), aligned('lowmq', mq=19),
                 aligned('duplicate', flag=1024), aligned('secondary', flag=256),
                 aligned('supplementary', flag=2048), aligned('qcfail', flag=512),
                 aligned('reverse', flag=16), aligned('unknown', base='N')]
        with pysam.AlignmentFile(str(self.bam('counts', reads)), 'rb') as bam:
            obs, reasons = collect(bam, 'I', 9, 20)
        result = support(obs, 'A', 'G', 13, True)
        self.assertEqual(result['depth'], 3)
        self.assertEqual(result['alt_count'], 2)
        self.assertEqual(result['other_count'], 1)
        self.assertEqual(result['alt_reverse'], 1)
        self.assertEqual(result['overlapping_mate_observations_collapsed'], 1)
        self.assertEqual(result['conflicting_eligible_fragments_excluded'], 1)
        self.assertEqual(reasons['excluded_alignment_flags'], 4)
        self.assertEqual(reasons['below_mapping_quality'], 1)
        self.assertEqual(support(obs, 'A','G',7,True)['depth'],4)
        self.assertEqual(support(obs,'A','G',13,False)['depth'],6)

    def test_deletion_refskip_and_softclip(self):
        for cigar in ('9M1D11M', '9M1N11M'):
            self.assertEqual(observation(aligned('gap', cigar=cigar),9)[1], 'deletion_or_refskip')
        read = aligned('clip', cigar='2S18M')
        obs, reason = observation(read,7)
        self.assertIsNone(reason)
        self.assertTrue(obs['soft_clipped'])
        self.assertEqual(obs['base'],'G')
        self.assertEqual(obs['aligned_end_distance'],7)
        self.assertIsNone(support([], 'A','G',13,True)['alt_fraction'])

    def test_end_to_end_and_reference_validation(self):
        left=self.bam('illumina',[aligned('pair',flag=65),aligned('pair',flag=129)])
        right=self.bam('ont',[aligned('read',quality=7)])
        out=self.root/'review'
        review(self.loci,self.ref,left,right,out)
        data=json.loads((out/'read_support.json').read_text())
        views=data['loci'][0]['views']
        self.assertEqual(views['illumina_bq13']['alt_count'],1)
        self.assertIsNone(views['ont_bq13']['alt_fraction'])
        self.assertEqual(views['ont_bq7']['alt_fraction'],1)
        self.assertEqual(data['visual_review_status'],'not_performed')
        self.assertEqual((out/'regions.bed').read_text(),'I\t9\t10\tshared:A>G\n')
        with (out/'read_support.tsv').open() as stream:
            self.assertEqual(len(list(csv.DictReader(stream,delimiter='\t'))),3)
        self.loci.write_text(self.loci.read_text().replace('\tA\tG\t','\tC\tG\t'))
        with self.assertRaisesRegex(ValueError,'Reference allele mismatch'):
            review(self.loci,self.ref,left,right,self.root/'invalid')
        self.assertFalse((self.root/'invalid'/'read_support.json').exists())

    def test_one_eligible_mate_and_empty_support(self):
        a=aligned('pair',base='G',quality=7,flag=65)
        b=aligned('pair',base='A',quality=20,flag=129)
        obs=[observation(r,9)[0] for r in (a,b)]
        self.assertEqual(support(obs,'A','G',13,True)['ref_count'],1)
        self.assertEqual(support(obs,'A','G',7,True)['depth'],0)
        self.assertEqual(support(obs,'A','G',7,True)['conflicting_eligible_fragments_excluded'],1)
        missing=aligned('missing')
        missing.query_qualities=None
        self.assertEqual(observation(missing,9)[1],'missing_sequence_or_quality')

    def test_locus_validation_and_portable_xml(self):
        self.loci.write_text(self.loci.read_text()+'I\t10\tA\tG\tshared\t0/1\t1/1\n')
        with self.assertRaisesRegex(ValueError,'Duplicate'):
            read_loci(self.loci)
        session=self.root/'session.xml'
        bam=self.root/'reads & mates.bam'
        igv_session(session,self.ref,{'Illumina & ONT':bam},'I:1-20')
        xml=ET.parse(session).getroot()
        self.assertEqual(xml.attrib['genome'],'reference.fa')
        r=xml.find('Resources/Resource')
        self.assertEqual(r.attrib['path'],'reads & mates.bam')
        self.assertEqual(r.attrib['index'],'reads & mates.bam.bai')

if __name__ == '__main__':
    unittest.main()
