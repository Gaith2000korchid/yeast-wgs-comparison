"""Native bcftools stage audit on real indexed synthetic BAMs."""
import importlib.util
from pathlib import Path
import random
import tempfile
import unittest

import pysam

spec = importlib.util.spec_from_file_location('caller_audit',
    Path(__file__).resolve().parents[1] / 'workflow/scripts/audit_regional_calls.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class CallerAuditTests(unittest.TestCase):
    def test_called_snp_is_removed_by_sample_depth_not_info_depth(self):
        with tempfile.TemporaryDirectory() as tmp:
            bundle = Path(tmp) / 'bundle'
            bundle.mkdir()
            rng = random.Random(12)
            sequence = ''.join(rng.choice('ACGT') for _ in range(600))
            reference = bundle / 'reference.fa'
            reference.write_text('>test\n' + sequence + '\n')
            pysam.faidx(str(reference))
            (bundle / 'extraction_windows.bed').write_text('test\t0\t600\n')
            targets = []
            for pos in (101, 301):
                ref = sequence[pos - 1]
                targets.append(('test', pos, ref, next(b for b in 'ACGT' if b != ref)))
            header = {'HD': {'SO': 'coordinate'}, 'SQ': [{'SN': 'test', 'LN': 600}],
                      'RG': [{'ID': 'x', 'SM': 'P11'}]}
            for tech in ('illumina', 'ont'):
                bam = bundle / f'P11.{tech}.regional.bam'
                with pysam.AlignmentFile(str(bam), 'wb', header=header) as stream:
                    for target, depth in zip(targets, (3, 12)):
                        _, pos, ref, alt = target
                        start = pos - 21
                        bases = list(sequence[start:start + 100])
                        bases[20] = alt
                        for n in range(depth):
                            read = pysam.AlignedSegment()
                            read.query_name = f'r{pos}_{n}'
                            read.query_sequence = ''.join(bases)
                            read.query_qualities = pysam.qualitystring_to_array('I' * 100)
                            read.reference_id = 0
                            read.reference_start = start
                            read.mapping_quality = 60
                            read.cigarstring = '100M'
                            read.set_tag('RG', 'x')
                            stream.write(read)
                pysam.index(str(bam))
            result = module.audit(bundle, Path(tmp) / 'audit', targets)
            for experiment in result['experiments'].values():
                records = experiment['records']
                self.assertEqual(len(records['raw']), 2)
                low = next(r for r in records['normalized'] if r['pos_1based'] == 101)
                self.assertGreater(low['qual'], 30)
                self.assertEqual(low['sample']['DP'], '3')
                self.assertEqual(low['rejection_reasons'], ['FMT/DP<8'])
                self.assertEqual([r['pos_1based'] for r in records['pass']], [301])
            with self.assertRaises(FileExistsError):
                module.audit(bundle, Path(tmp) / 'audit', targets)


if __name__ == '__main__':
    unittest.main()
