"""Synthetic cache contracts and direct consumption by the original anchor indexer."""
import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
from common import COUNTS, QCOLS, read, sha, write
from run import execute

GUARDS = dict(reserve_memory_gib=0, reserve_disk_gib=0, worker_gib=8)


def row(chrom=1, ps='001', pos=100):
    return dict(sample='sample01', tissue='Heart', chrom=str(chrom), cpg_pos1=str(pos), ps=ps,
        genotype_disrupted='0', **dict(zip(COUNTS, ['5', '5', '5', '5'])),
        H1_depth='10', H2_depth='10', H1_methylation='0.5', H2_methylation='0.5', delta_M='0',
        test_status='TESTED', fisher_two_sided_p='1', **{q: '1' for q in QCOLS},
        annotation_gene_id='synthetic_gene', annotation_primary_feature='PROMOTER')


def fixture(root, records=None):
    records = records or [row(c) for c in range(1, 19)] + [row(ps='1')]
    p = root / 'annotated.tsv.gz'
    with gzip.open(p, 'wt', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(records[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(records)
    spec = dict(sample='sample01', tissue='Heart', path=p.name, sha256=sha(p), rows=len(records))
    config = root / 'inputs.json'
    write(config, dict(units=[dict(sample='sample01', tissue='Heart', annotated=spec)]))
    return config, p


class CacheTests(unittest.TestCase):
    def test_preservation_restart_and_native_consumer(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); config, source = fixture(root)
            registry = execute(config, root / 'cache', GUARDS)
            self.assertEqual(registry, execute(config, root / 'cache', GUARDS))
            unit = registry['units'][0]; cp = Path(unit['cache']['path'])
            self.assertEqual(read(cp / 'phase_dictionary.json'), ['001', '1'])
            a = np.load(cp / 'chr1.npy', allow_pickle=False)
            self.assertEqual(a['source_row'].tolist(), [1, 19])
            self.assertEqual(a['p'].tolist(), [1, 1])  # Non-significant rows retained.
            for name in ['a', 'b', 'c', 'd']:
                self.assertEqual(a[name].tolist(), [5, 5])
            for name in ['q_unit', 'q_tissue', 'q_all']:
                self.assertEqual(a[name].tolist(), [1, 1])
            with gzip.open(cp / 'chr1.source.tsv.gz', 'rt', encoding='utf-8') as f:
                saved = list(csv.DictReader(f, delimiter='\t'))
            self.assertEqual(saved[1]['ps'], '1')
            self.assertEqual(saved[0]['annotation_gene_id'], 'synthetic_gene')
            unit['raw'] = dict(unit['annotated'])
            job = dict(spec=unit, parent=str(root / 'availability'), identity='synthetic-integration', name='sample01-Heart')
            write(root / 'job.json', job)
            package = Path(__file__).resolve().parents[1] / 'asm_allele_prep'
            # Only resource admission is disabled for this tiny synthetic run.
            code = ('import sys,json;sys.path.insert(0,sys.argv[1]);import indexing;'
                    'indexing.guard=lambda *a,**k:None;'
                    'print(indexing.unit_job(json.load(open(sys.argv[2],encoding="utf-8"))))')
            r = subprocess.run([sys.executable, '-X', 'utf8', '-c', code, str(package), str(root / 'job.json')],
                               capture_output=True, text=True, encoding='utf-8')
            self.assertEqual(r.returncode, 0, r.stderr)
            out = Path(r.stdout.strip().splitlines()[-1])
            self.assertEqual(read(out / 'summary.json')['tested_rows'], 19)

    def test_duplicate_exact_ps_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); config, _ = fixture(root, [row(), row()])
            with self.assertRaisesRegex(ValueError, 'Duplicate same-CpG'):
                execute(config, root / 'cache', GUARDS)

    def test_count_arithmetic_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); rec = row(); rec['H1_depth'] = '11'
            config, _ = fixture(root, [rec])
            with self.assertRaisesRegex(ValueError, 'Depth differs'):
                execute(config, root / 'cache', GUARDS)

    def test_source_hash_change_rejected_on_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); config, source = fixture(root)
            execute(config, root / 'cache', GUARDS)
            source.write_bytes(source.read_bytes() + b'changed')
            with self.assertRaisesRegex(ValueError, 'Input identity changed'):
                execute(config, root / 'cache', GUARDS)

    def test_relative_input_paths_from_other_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); config, _ = fixture(root)
            elsewhere = root / 'elsewhere'; elsewhere.mkdir()
            command = [sys.executable, '-X', 'utf8', str(Path(__file__).with_name('run.py')),
                '--config', str(config), '--output', str(root / 'cache'),
                '--reserve-memory-gib', '0', '--reserve-disk-gib', '0', '--worker-memory-gib', '8']
            r = subprocess.run(command, cwd=elsewhere, capture_output=True, text=True, encoding='utf-8')
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(read(root / 'cache/validation.json')['rows'], 19)


if __name__ == '__main__':
    unittest.main()
