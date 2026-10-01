"""Check paired sampling and historical ONT molecular independence."""
import gzip
import io
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
from prepare_real_reads import records, subsample_pairs, combine_ont
from download_reads import verify
import hashlib


def fastq(name, seq='ACGT'):
    return f'@{name}\n{seq}\n+\n' + 'I' * len(seq) + '\n'


class PreparationTests(unittest.TestCase):
    def test_sampling_retains_pairs_and_is_byte_reproducible(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            for mate in [1, 2]:
                with gzip.open(p / f'in{mate}.gz', 'wt') as f:
                    f.write(''.join(fastq(f'read{i}/{mate}') for i in range(200)))
            first = subsample_pairs(p/'in1.gz', p/'in2.gz', p/'a.gz', p/'b.gz', .25, 42)
            second = subsample_pairs(p/'in1.gz', p/'in2.gz', p/'c.gz', p/'d.gz', .25, 42)
            self.assertEqual(first, second)
            self.assertEqual(first['input_pairs'], 200)
            self.assertTrue(0 < first['selected_pairs'] < 200)
            self.assertEqual((p/'a.gz').read_bytes(), (p/'c.gz').read_bytes())
            with gzip.open(p/'a.gz') as a, gzip.open(p/'b.gz') as b:
                left=list(records(a));right=list(records(b))
            self.assertEqual(len(left), first['selected_pairs'])
            self.assertEqual([r[0].replace(b'/1',b'') for r in left], [r[0].replace(b'/2',b'') for r in right])

    def test_mismatched_mates_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            for mate,name in [(1,'one'),(2,'two')]:
                with gzip.open(p/f'{mate}.gz','wt') as f:f.write(fastq(name))
            with self.assertRaisesRegex(ValueError,'synchronized'):
                subsample_pairs(p/'1.gz',p/'2.gz',p/'a.gz',p/'b.gz',1,1)

    def test_ont_retains_one_consensus_and_reports_read_classes(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            text=''.join(fastq(f'ERR1.{i} ch1_file1_{kind}_fail_runA/1') for i,kind in enumerate(['template','complement','twodirections']))
            text+=fastq('ERR1.4 ch1_file2_twodirections_pass_runA/1')
            with gzip.open(p/'in.gz','wt') as f:f.write(text)
            result=combine_ont([p/'in.gz'],p/'out.gz')[0]
            self.assertEqual(result['archive_reads'],4)
            self.assertEqual(result['selected_2d_reads'],2)
            with gzip.open(p/'out.gz') as f:
                self.assertTrue(all(b'twodirections' in r[0] for r in records(f)))
            with self.assertRaisesRegex(ValueError,'Duplicate'):
                combine_ont([p/'in.gz',p/'in.gz'],p/'bad.gz')

    def test_truncated_fastq_and_incorrect_checksum_rejected(self):
        with self.assertRaises(ValueError):list(records(io.BytesIO(b'@x\nACGT\n+\n')))
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x';p.write_bytes(b'hello')
            self.assertTrue(verify(p,{'bytes':5,'md5':hashlib.md5(b'hello').hexdigest()}))
            self.assertFalse(verify(p,{'bytes':5,'md5':'0'*32}))


if __name__=='__main__':unittest.main()
