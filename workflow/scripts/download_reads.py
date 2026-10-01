"""Download manifest FASTQs with atomic writes and archive byte/MD5 validation."""
import argparse
import hashlib
import json
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def verify(path, entry):
    if not path.is_file() or path.stat().st_size != entry["bytes"]:
        return False
    digest = hashlib.md5()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest() == entry["md5"]


def download(entry, out):
    path = Path(out) / entry["filename"]
    if path.name != entry["filename"] or entry["filename"] in {".", ".."}:
        raise ValueError("Manifest filenames must be basenames")
    if not entry["url"].startswith("https://ftp.sra.ebi.ac.uk/"):
        raise ValueError("Expected an ENA HTTPS archive URL")
    if verify(path, entry):
        print(f"Already verified: {path.name}", flush=True)
        return
    if path.exists():
        raise ValueError(f"Existing file failed checksum/size validation: {path}")
    partial = path.with_name(path.name + ".part")
    for attempt in range(1, 4):
        digest = hashlib.md5()
        total = next_report = 0
        try:
            print(f"Downloading {path.name}: {entry['bytes'] / 1e6:.1f} MB (attempt {attempt})", flush=True)
            request = urllib.request.Request(entry["url"], headers={"User-Agent": "yeast-wgs-comparison/0.2"})
            with urllib.request.urlopen(request, timeout=60) as response, partial.open("wb") as stream:
                while block := response.read(1024 * 1024):
                    stream.write(block)
                    digest.update(block)
                    total += len(block)
                    if total >= next_report:
                        print(f"{path.name}: {total / 1e6:.1f} / {entry['bytes'] / 1e6:.1f} MB", flush=True)
                        next_report = total + 100 * 1024 * 1024
            if total != entry["bytes"] or digest.hexdigest() != entry["md5"]:
                raise ValueError(f"Archive byte count/MD5 mismatch for {path.name}")
            partial.replace(path)
            print(f"Verified: {path.name} ({digest.hexdigest()})", flush=True)
            return
        except Exception as error:
            print(f"Attempt {attempt} failed for {path.name}: {error}", flush=True)
            if attempt == 3:
                raise
            time.sleep(2 * attempt)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", required=True)
    p.add_argument("--out", default="data/real/raw")
    p.add_argument("--workers", type=int, default=2)
    a = p.parse_args()
    if not 1 <= a.workers <= 4:
        p.error("Use 1 to 4 workers")
    entries = json.loads(Path(a.manifest).read_text())["files"]
    Path(a.out).mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=a.workers) as executor:
        list(executor.map(lambda entry: download(entry, a.out), entries))
