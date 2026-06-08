#!/usr/bin/env python3
"""
Pobiera dane multi-omics z BigQuery używając ARRAY_AGG (jeden wiersz = jeden pacjent).
"""

import pandas as pd
from google.cloud import bigquery
import os

PROJECT_ID = "isb-cgc-bq"

PROJECTS = {
    "TCGA": ["cnv", "rna", "mirna", "methylation"],
    "CPTAC": ["cnv", "rna", "mirna"],
    "TARGET": ["cnv", "rna", "mirna"],
    "CDDP_EAGLE": ["cnv", "rna"],
    "CGCI": ["cnv", "rna"],
    "HCMI": ["cnv", "rna"],
    "MP2PRT": ["cnv", "rna"],
    "REBC": ["cnv", "rna"],
    "APOLLO": ["rna"],
    "BEATAML1_0": ["rna"],
    "CMI": ["rna"],
    "CTSP": ["rna"],
    "EXC_RESPONDERS": ["rna"],
    "MMRF": ["rna"],
    "NCICCR": ["rna"],
    "OHSU": ["rna"],
    "ORGANOID": ["rna"],
    "WCDT": ["rna"],
}

def get_table_name(project, modality):
    tables = {
        "cnv":        f"isb-cgc-bq.{project}.copy_number_gene_level_hg38_gdc_current",
        "rna":        f"isb-cgc-bq.{project}.RNAseq_hg38_gdc_current",
        "mirna":      f"isb-cgc-bq.{project}.miRNAseq_hg38_gdc_current",
        "methylation":f"isb-cgc-bq.{project}.DNA_methylation_hg38_gdc_current",
    }
    return tables[modality]

def load_gene_lists():
    df = pd.read_csv('https://raw.githubusercontent.com/marcin119a/data/refs/heads/main/dataset_columns.csv')
    return {m: df[df['modality'] == m]['column_name'].unique().tolist()
            for m in ['cnv', 'rna', 'mirna', 'methylation']}

def fetch_cnv(client, table, gene_list):
    genes_sql = "', '".join(gene_list)
    query = f"""
    SELECT
        case_barcode,
        ARRAY_AGG(gene_name   IGNORE NULLS ORDER BY gene_name) AS labels,
        ARRAY_AGG(copy_number IGNORE NULLS ORDER BY gene_name) AS values
    FROM `{table}`
    WHERE gene_name IN ('{genes_sql}')
    GROUP BY case_barcode
    """
    try:
        return client.query(query).to_dataframe()
    except Exception as e:
        print(f"    ⚠ Błąd CNV: {e}")
        return None

def fetch_rna(client, table, gene_list):
    genes_sql = "', '".join(gene_list)
    query = f"""
    SELECT
        case_barcode,
        ARRAY_AGG(gene_name        IGNORE NULLS ORDER BY gene_name) AS labels,
        ARRAY_AGG(fpkm_unstranded  IGNORE NULLS ORDER BY gene_name) AS values
    FROM `{table}`
    WHERE gene_name IN ('{genes_sql}')
    GROUP BY case_barcode
    """
    try:
        return client.query(query).to_dataframe()
    except Exception as e:
        print(f"    ⚠ Błąd RNA: {e}")
        return None

def fetch_mirna(client, table, mirna_list):
    mirnas_sql = "', '".join(mirna_list)
    query = f"""
    SELECT
        case_barcode,
        ARRAY_AGG(mirna_id                        IGNORE NULLS ORDER BY mirna_id) AS labels,
        ARRAY_AGG(reads_per_million_miRNA_mapped  IGNORE NULLS ORDER BY mirna_id) AS values
    FROM `{table}`
    WHERE mirna_id IN ('{mirnas_sql}')
    GROUP BY case_barcode
    """
    try:
        return client.query(query).to_dataframe()
    except Exception as e:
        print(f"    ⚠ Błąd miRNA: {e}")
        return None

def fetch_methylation(client, table, probe_list):
    probes_sql = "', '".join(probe_list)
    query = f"""
    SELECT
        case_barcode,
        ARRAY_AGG(probe_id    IGNORE NULLS ORDER BY probe_id) AS labels,
        ARRAY_AGG(beta_value  IGNORE NULLS ORDER BY probe_id) AS values
    FROM `{table}`
    WHERE probe_id IN ('{probes_sql}')
    GROUP BY case_barcode
    """
    try:
        return client.query(query).to_dataframe()
    except Exception as e:
        print(f"    ⚠ Błąd metylacji: {e}")
        return None

def main():
    genes = load_gene_lists()
    print(f"Wczytano: {sum(len(v) for v in genes.values())} genów/sond\n")

    client = bigquery.Client()
    os.makedirs("data_hg38", exist_ok=True)

    fetchers = {
        'cnv':        (fetch_cnv,        'cnv'),
        'rna':        (fetch_rna,        'rna'),
        'mirna':      (fetch_mirna,      'mirna'),
        'methylation':(fetch_methylation,'methylation'),
    }

    for project, modalities in PROJECTS.items():
        print(f"\n{'='*60}")
        print(f"Projekt: {project} {modalities}")
        print('='*60)

        for modality in modalities:
            print(f"  Pobieranie {modality.upper()}...")
            fetch_fn, gene_key = fetchers[modality]
            table = get_table_name(project, modality)
            df = fetch_fn(client, table, genes[gene_key])

            if df is not None and not df.empty:
                print(f"    ✓ {modality.upper()}: {df.shape}")
                path = f"data_hg38/{project}_{modality}_hg38.parquet"
                df.to_parquet(path, index=False)
                print(f"  ✓ Zapisano: {path}")
            else:
                print(f"    ✗ {modality.upper()}: brak danych")

    print(f"\n{'='*60}")
    print("✓ Gotowe! Dane zapisane w folderze data_hg38/")
    print('='*60)

if __name__ == "__main__":
    main()
