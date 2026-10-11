from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: str | Path, rows: list[dict], columns: list[str] | None = None) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    columns = columns or (list(rows[0]) if rows else [])
    if not columns:
        raise ValueError("Output columns are required for an empty table")
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return target


def read_json(path: str | Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def nonmissing(value: object) -> bool:
    return str(value if value is not None else "").strip().lower() not in {
        "", "na", "nan", "none", "null", "-", "not tested"
    }


def float_or_none(value: object) -> float | None:
    if not nonmissing(value):
        return None
    try:
        result = float(str(value))
    except (TypeError, ValueError):
        return None
    import math
    return result if math.isfinite(result) else None
