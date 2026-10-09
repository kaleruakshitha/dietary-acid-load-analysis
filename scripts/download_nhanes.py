"""Download the verified core NHANES XPT files.

Run from the repository root:

    python scripts/download_nhanes.py

Files are written only when absent unless ``--overwrite`` is supplied. A
SHA-256 manifest is generated so later pipeline runs can verify source bytes.
"""

from __future__ import annotations

import argparse
import hashlib
import time
import urllib.request
from pathlib import Path


BASE_URL = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles"
FILES = [
    "DEMO_L.xpt",
    "DR1TOT_L.xpt",
    "DR2TOT_L.xpt",
    "BIOPRO_L.xpt",
    "BMX_L.xpt",
    "ALB_CR_L.xpt",
    "SMQ_L.xpt",
    "DIQ_L.xpt",
    "KIQ_U_L.xpt",
    "DR1IFF_L.xpt",
    "DR2IFF_L.xpt",
    "DRXFCD_L.xpt",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    destination = Path("data/raw/nhanes")
    destination.mkdir(parents=True, exist_ok=True)

    manifest_rows = ["file_name,sha256"]
    for file_name in FILES:
        target = destination / file_name
        if args.overwrite or not target.exists():
            request = urllib.request.Request(
                f"{BASE_URL}/{file_name}", headers={"User-Agent": "Mozilla/5.0"}
            )
            last_error = None
            for attempt in range(4):
                try:
                    with urllib.request.urlopen(request, timeout=120) as response:
                        target.write_bytes(response.read())
                    last_error = None
                    break
                except Exception as error:  # network errors vary by platform
                    last_error = error
                    if attempt < 3:
                        time.sleep(2**attempt)
            if last_error is not None:
                raise last_error
        manifest_rows.append(f"{file_name},{sha256(target)}")

    (destination / "sha256_manifest.csv").write_text(
        "\n".join(manifest_rows) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()