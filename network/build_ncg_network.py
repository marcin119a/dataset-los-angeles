"""
Build the NCG Network from DeepGraphMut (DGM) paper.

Steps:
  1. Collect gene symbols from NCG cancer-driver + healthy-driver TSVs.
  2. Parse PCNet CX file to get node id -> gene symbol and all edges.
  3. Keep only nodes whose gene symbol appears in the NCG gene set.
  4. Keep only edges where both endpoints are retained nodes.
  5. Write edges as a two-column TSV (gene_a, gene_b).
"""

import json
import csv
import sys
from pathlib import Path

DATA_DIR = Path(__file__).parent

CANCER_DRIVERS = DATA_DIR / "NCG_cancerdrivers_annotation_supporting_evidence.tsv"
HEALTHY_DRIVERS = DATA_DIR / "NCG_healthydrivers_annotation_supporting_evidence.tsv"
PCNET_CX = DATA_DIR / "Parsimonious Composite Network (PCNet).cx"
OUTPUT_EDGES = DATA_DIR / "NCG_network_edges.tsv"
OUTPUT_NODES = DATA_DIR / "NCG_network_nodes.txt"


def load_ncg_genes(*tsv_paths: Path) -> set[str]:
    genes: set[str] = set()
    for path in tsv_paths:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                symbol = row["symbol"].strip()
                if symbol:
                    genes.add(symbol)
    return genes


def parse_pcnet(cx_path: Path) -> tuple[dict[int, str], list[tuple[int, int]]]:
    with open(cx_path) as fh:
        cx = json.load(fh)

    id_to_gene: dict[int, str] = {}
    raw_edges: list[tuple[int, int]] = []

    for section in cx:
        if "nodes" in section:
            for node in section["nodes"]:
                id_to_gene[node["@id"]] = node["n"]
        if "edges" in section:
            for edge in section["edges"]:
                raw_edges.append((edge["s"], edge["t"]))

    return id_to_gene, raw_edges


def build_network(
    ncg_genes: set[str],
    id_to_gene: dict[int, str],
    raw_edges: list[tuple[int, int]],
) -> tuple[set[str], list[tuple[str, str]]]:
    # Node IDs whose gene symbol is in the NCG gene set
    valid_ids = {nid for nid, sym in id_to_gene.items() if sym in ncg_genes}

    edges: list[tuple[str, str]] = []
    seen: set[frozenset[str]] = set()
    for s, t in raw_edges:
        if s in valid_ids and t in valid_ids:
            gene_s = id_to_gene[s]
            gene_t = id_to_gene[t]
            key = frozenset({gene_s, gene_t})
            if key not in seen and gene_s != gene_t:
                seen.add(key)
                edges.append((gene_s, gene_t))

    retained_genes = {g for pair in edges for g in pair}
    return retained_genes, edges


def main() -> None:
    print("Loading NCG gene sets …")
    ncg_genes = load_ncg_genes(CANCER_DRIVERS, HEALTHY_DRIVERS)
    print(f"  NCG genes (unique symbols): {len(ncg_genes)}")

    print("Parsing PCNet CX …")
    id_to_gene, raw_edges = parse_pcnet(PCNET_CX)
    print(f"  PCNet nodes: {len(id_to_gene)}, edges: {len(raw_edges)}")

    print("Filtering to NCG subnetwork …")
    retained_genes, edges = build_network(ncg_genes, id_to_gene, raw_edges)
    print(f"  Retained nodes: {len(retained_genes)}, edges: {len(edges)}")

    print(f"Writing edges → {OUTPUT_EDGES}")
    with open(OUTPUT_EDGES, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["gene_a", "gene_b"])
        writer.writerows(edges)

    print(f"Writing node list → {OUTPUT_NODES}")
    OUTPUT_NODES.write_text("\n".join(sorted(retained_genes)) + "\n")

    print("Done.")
    print(f"  Nodes: {len(retained_genes)}  (paper reports 3217)")
    print(f"  Edges: {len(edges)}")


if __name__ == "__main__":
    main()
