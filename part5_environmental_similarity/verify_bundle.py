"""Verify delivered file integrity using only the Python standard library."""
from pathlib import Path
import csv
import hashlib
root = Path(__file__).resolve().parent
with (root / "bundle_checksums.csv").open(newline="") as stream:
    rows = list(csv.DictReader(stream))
for row in rows:
    path = root / row["path"]
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
        raise SystemExit("Changed or missing: " + row["path"])
print(f"Verified {len(rows)} delivered files.")
