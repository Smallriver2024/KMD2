from __future__ import annotations

from pathlib import Path

from .io import float_or_none


def write_gct(path: str | Path, rows: list[dict], sample_ids: list[str]) -> Path:
    if not rows:
        raise ValueError("Refusing to write an empty PTM-SEA GCT")
    if len(sample_ids) != len(set(sample_ids)) or not sample_ids:
        raise ValueError("GCT sample IDs must be nonempty and unique")
    names = [r["Name"] for r in rows]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate PTM-SEA row identifiers")
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="") as handle:
        handle.write("#1.2\n")
        handle.write(f"{len(rows)}\t{len(sample_ids)}\n")
        handle.write("\t".join(["Name", "Description", *sample_ids]) + "\n")
        for row in rows:
            vals = ["NA" if row.get(s) is None else f"{float(row[s]):.6f}" for s in sample_ids]
            handle.write("\t".join([row["Name"], row["Description"], *vals]) + "\n")
    return target


def read_gct(path: str | Path) -> tuple[list[str], list[dict]]:
    with Path(path).open("r", encoding="utf-8-sig") as handle:
        version = handle.readline().strip()
        if version not in {"#1.2", "#1.3"}:
            raise ValueError(f"Unsupported GCT version: {version}")
        dimensions = handle.readline().strip().split("\t")
        header = handle.readline().rstrip("\n\r").split("\t")
        if version == "#1.2":
            sample_ids = header[2:]
        else:
            if len(dimensions) != 4:
                raise ValueError("GCT 1.3 requires row/column and metadata dimensions")
            n_samples, n_row_metadata, n_column_metadata = map(int, dimensions[1:])
            if len(header) != 1 + n_row_metadata + n_samples:
                raise ValueError("GCT 1.3 header width differs from declared dimensions")
            sample_ids = header[1 + n_row_metadata:]
            for _ in range(n_column_metadata):
                metadata = handle.readline().rstrip("\n\r").split("\t")
                if len(metadata) != len(header):
                    raise ValueError("GCT column metadata width differs from header")
        if not sample_ids or len(sample_ids) != len(set(sample_ids)):
            raise ValueError("GCT sample headers invalid")
        result = []
        for line in handle:
            cells = line.rstrip("\n\r").split("\t")
            if len(cells) != len(header):
                raise ValueError("GCT row width differs from header")
            row = {"signature": cells[0]}
            for sample in sample_ids:
                row[sample] = float_or_none(cells[header.index(sample)])
            result.append(row)
    if int(dimensions[0]) != len(result):
        raise ValueError("GCT declared row count differs from content")
    if version == "#1.2" and int(dimensions[1]) != len(sample_ids):
        raise ValueError("GCT declared sample count differs from content")
    return sample_ids, result
