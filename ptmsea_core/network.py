from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from .io import nonmissing, read_tsv, sha256
from .config import normalize_symbols, target_aliases


def prepare_ppi_sources(config_rows: list[dict], repo: Path) -> list[dict]:
    """Consume licensed, locally obtained source tables with explicit column maps."""
    edges = []
    for source in config_rows:
        if source["active"].lower() != "true":
            continue
        path = Path(source["source_file"])
        if not path.is_absolute():
            path = repo / path
        if not path.exists():
            raise FileNotFoundError(f"{source['source_name']}: {path}")
        if source.get("sha256") and sha256(path).lower() != source["sha256"].lower():
            raise ValueError(f"PPI source hash mismatch: {source['source_name']}")
        columns = read_tsv(path)
        for row in columns:
            a, b = row[source["protein_a_col"]].strip(), row[source["protein_b_col"]].strip()
            if a and b:
                edges.append({
                    "protein_a": a, "protein_b": b,
                    "protein_a_id_type": source["id_type"], "protein_b_id_type": source["id_type"],
                    "source_database": source["source_name"], "source_version": source["source_version"],
                    "interaction_category": source["edge_type"],
                    "interaction_type": row.get(source.get("interaction_type_col", ""), ""),
                    "pmid": row.get(source.get("pmid_col", ""), ""),
                })
    return edges


def build_reference_network(edges: list[dict], id_to_gene: dict[str, str] | None = None) -> list[dict]:
    """Aggregate one-hop undirected gene edges and preserve all source labels."""
    mapping = id_to_gene or {}
    if any(edge["protein_a_id_type"] != "gene_symbol" or edge["protein_b_id_type"] != "gene_symbol"
           for edge in edges) and not mapping:
        raise ValueError("Non-gene PPI identifiers require an explicit approved-symbol map")
    pairs: dict[tuple[str, str], dict[str, set[str] | int]] = {}
    for edge in edges:
        a = mapping.get(edge["protein_a"], edge["protein_a"] if edge["protein_a_id_type"] == "gene_symbol" else "")
        b = mapping.get(edge["protein_b"], edge["protein_b"] if edge["protein_b_id_type"] == "gene_symbol" else "")
        a, b = a.strip().upper(), b.strip().upper()
        if not a or not b or a == b:
            continue
        key = tuple(sorted((a, b)))
        item = pairs.setdefault(key, {"sources": set(), "categories": set(), "types": set(), "pmids": set(), "n_edges": 0})
        item["sources"].add(edge["source_database"])
        item["categories"].add(edge.get("interaction_category", ""))
        if edge.get("interaction_type"):
            item["types"].add(edge["interaction_type"])
        if edge.get("pmid"):
            item["pmids"].update(value.strip() for value in str(edge["pmid"]).replace("|", ";").split(";") if nonmissing(value))
        item["n_edges"] += 1
    return [
        {"gene_a": a, "gene_b": b, "source_database": ";".join(sorted(v["sources"])),
         "interaction_category": ";".join(sorted(v["categories"])),
         "interaction_type": ";".join(sorted(v["types"])),
         "pmid": ";".join(sorted(x for x in v["pmids"] if x)),
         "n_source_edges": v["n_edges"]}
        for (a, b), v in sorted(pairs.items())
    ]


def annotate_network_support(evidence_rows: list[dict], network_rows: list[dict],
                             aliases: dict[str, list[str]] | None = None) -> list[dict]:
    """Keep the record universe and favor strict drug-target support over kinase support."""
    network = {frozenset((r["gene_a"].upper(), r["gene_b"].upper())): r for r in network_rows}
    aliases = target_aliases() if aliases is None else aliases
    out = []
    for evidence in evidence_rows:
        genes = normalize_symbols(evidence.get("regulated_gene", ""), aliases)
        gene = genes[0] if len(genes) == 1 else ""
        targets = normalize_symbols(evidence.get("drug_target", ""), aliases, require_approved=True)
        centers = normalize_symbols(evidence.get("center_kinase", ""), aliases)
        result = dict(evidence)
        result.update(regulated_gene_approved=gene, drug_target_approved=";".join(targets))
        result.update(network_support_status="none", network_match_type="",
                      network_pair="", network_sources="", network_pmids="",
                      network_match_count=0)
        if gene and gene in targets:
            result.update(network_support_status="direct_target", network_match_type="candidate_equals_drug_target")
        else:
            for role, symbols in (("drug_target", targets), ("center_kinase", centers)):
                matches = [(symbol, network[frozenset((gene, symbol))]) for symbol in symbols
                           if gene and gene != symbol and frozenset((gene, symbol)) in network]
                if matches:
                    symbol, edge = matches[0]
                    result.update(network_support_status="direct_partner",
                                  network_match_type=f"one_hop_{role}",
                                  network_pair="--".join(sorted((gene, symbol))),
                                  network_sources=edge["source_database"],
                                  network_pmids=edge.get("pmid", ""),
                                  network_match_count=int(edge["n_source_edges"]))
                    break
        out.append(result)
    if len(out) != len(evidence_rows):
        raise AssertionError("Network annotation changed the phosphoproteomics record count")
    return out
