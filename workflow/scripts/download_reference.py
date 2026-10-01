"""Fetch the pinned Ensembl reference and check both recorded SHA256 digests."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import urllib.request


def digest(path):
    value = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def fetch(manifest, output):
    metadata = json.loads(Path(manifest).read_text())
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.is_file() and digest(output) == metadata['fasta_sha256']:
        print('Reference already verified')
        return
    archive = output.with_name(output.name + '.gz.part')
    fasta = output.with_name(output.name + '.part')
    with urllib.request.urlopen(metadata['url'], timeout=60) as source, archive.open('wb') as target:
        shutil.copyfileobj(source, target)
    if digest(archive) != metadata['archive_sha256']:
        raise ValueError('Reference archive SHA256 mismatch')
    with gzip.open(archive, 'rb') as source, fasta.open('wb') as target:
        shutil.copyfileobj(source, target)
    if digest(fasta) != metadata['fasta_sha256']:
        raise ValueError('Uncompressed reference SHA256 mismatch')
    fasta.replace(output)
    archive.unlink()
    print('Pinned reference downloaded and verified')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', default='config/P11.reference.json')
    p.add_argument('--out', default='data/real/reference.fa')
    a = p.parse_args()
    fetch(a.manifest, a.out)
