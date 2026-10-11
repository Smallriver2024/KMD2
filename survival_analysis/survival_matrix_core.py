"""Core helpers for summarising KinaseMD survival-analysis outputs.

Adapted from the survival-result summarisation workflow.
The original log-rank p-values are read from the upstream OncoLnc-derived
result file; this module intentionally contains no plotting or source data.
"""

from __future__ import annotations

import csv
from collections import OrderedDict
from collections.abc import Iterable
from pathlib import Path


def split_sql_tuple(tuple_text: str) -> list[str | None]:
    """Split one SQL VALUES tuple while respecting quoted commas and escapes."""
    fields: list[str | None] = []
    current: list[str] = []
    in_quote = escape = quoted_field = False
    for char in tuple_text:
        if escape:
            current.append(char)
            escape = False
        elif in_quote and char == "\\":
            escape = True
        elif char == "'":
            in_quote = not in_quote
            quoted_field = True
        elif char == "," and not in_quote:
            value = "".join(current).strip()
            fields.append(value if quoted_field or value.upper() != "NULL" else None)
            current, quoted_field = [], False
        else:
            current.append(char)
    value = "".join(current).strip()
    fields.append(value if quoted_field or value.upper() != "NULL" else None)
    return fields


def extract_sql_tuples(text: str) -> list[str]:
    """Extract top-level SQL VALUES tuples without splitting inside quotes."""
    tuples: list[str] = []
    current: list[str] = []
    in_quote = escape = False
    depth = 0
    for char in text:
        if escape:
            if depth:
                current.append(char)
            escape = False
        elif in_quote and char == "\\":
            if depth:
                current.append(char)
            escape = True
        elif char == "'":
            if depth:
                current.append(char)
            in_quote = not in_quote
        elif not in_quote and char == "(":
            if depth:
                current.append(char)
            else:
                current = []
            depth += 1
        elif not in_quote and char == ")" and depth:
            depth -= 1
            if depth:
                current.append(char)
            else:
                tuples.append("".join(current))
        elif depth:
            current.append(char)
    return tuples


def gene_to_entrez_from_sql(sql_text: str) -> dict[str, str]:
    """Map kinase symbols to Entrez identifiers from a kinase SQL export."""
    mapping: dict[str, str] = {}
    for tuple_text in extract_sql_tuples(sql_text):
        fields = split_sql_tuple(tuple_text)
        if len(fields) >= 7 and fields[1] and fields[6]:
            mapping[str(fields[1])] = str(fields[6])
    return mapping


def parse_logrank_records(lines: Iterable[str]) -> tuple[dict[tuple[str, str], float], list[str]]:
    """Parse ``GENE_CANCER\tLogrank p-value=<p>`` records."""
    pvalues: dict[tuple[str, str], float] = {}
    cancers: OrderedDict[str, None] = OrderedDict()
    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            continue
        key, p_text = line.split("\t", 1)
        gene, cancer = key.rsplit("_", 1)
        pvalues[(gene, cancer)] = float(p_text.removeprefix("Logrank p-value="))
        cancers[cancer] = None
    return pvalues, list(cancers)


def read_group_counts(path: Path) -> tuple[int, int, int]:
    """Return total, low-expression, and high-expression sample counts."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        expected = ["Patient", "Days", "Status", "Expression", "Group"]
        if reader.fieldnames != expected:
            raise ValueError(f"Unexpected CSV header in {path}: {reader.fieldnames}")
        rows = list(reader)
    return len(rows), sum(row["Group"] == "Low" for row in rows), sum(row["Group"] == "High" for row in rows)


def make_result_row(gene: str, cancer: str, p_value: float, counts: tuple[int, int, int]) -> dict[str, object]:
    """Build the concise kinase-by-cancer survival result record."""
    total, low, high = counts
    significance = "Sig" if p_value < 0.05 else "NS"
    return {
        "Kinase": gene,
        "Cancer": cancer,
        "p_value": p_value,
        "significance": significance,
        "sample_total": total,
        "sample_low": low,
        "sample_high": high,
        "result": f"{significance}; p={p_value}; n={total}; Low={low}; High={high}",
    }
