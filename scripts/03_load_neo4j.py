"""
Charge les CSVs traités (data/processed/) dans Neo4j :
  - (:Protein {name})
  - (:Protein)-[:INTERACTS_WITH {score}]->(:Protein)
  - (:Drug {chembl_id, name, warning_class, warning_year})
  - (:Drug)-[:TARGETS {action_type}]->(:Protein)

Idempotent : utilise MERGE partout, donc on peut relancer sans dupliquer.
"""
import os
from pathlib import Path

import pandas as pd
from dotenv import load_dotenv
from neo4j import GraphDatabase
from tqdm import tqdm

load_dotenv()

URI = os.getenv("NEO4J_URI", "bolt://localhost:7687")
USER = os.getenv("NEO4J_USER", "neo4j")
PASSWORD = os.getenv("NEO4J_PASSWORD", "toxkg-dev-pw")

PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"

BATCH_SIZE = 2000


def batches(df: pd.DataFrame, size: int):
    for i in range(0, len(df), size):
        yield df.iloc[i : i + size]


def load_ppi(session):
    df = pd.read_csv(PROC_DIR / "ppi_edges.csv")
    print(f"[ppi] {len(df)} interactions à charger")

    query = """
    UNWIND $rows AS row
    MERGE (a:Protein {name: row.source})
    MERGE (b:Protein {name: row.target})
    MERGE (a)-[r:INTERACTS_WITH]-(b)
    SET r.score = row.score
    """
    for batch in tqdm(list(batches(df, BATCH_SIZE)), desc="loading PPI"):
        session.run(query, rows=batch.to_dict("records"))


def load_drugs_and_targets(session):
    drugs = pd.read_csv(PROC_DIR / "toxic_drugs.csv")
    targets = pd.read_csv(PROC_DIR / "drug_targets.csv")

    print(f"[drugs] {len(drugs)} médicaments toxiques")
    query_drugs = """
    UNWIND $rows AS row
    MERGE (d:Drug {chembl_id: row.chembl_id})
    SET d.name = row.pref_name,
        d.toxicity_reason = row.toxicity_reason,
        d.warning_class = row.warning_class,
        d.warning_country = row.warning_country,
        d.warning_year = row.warning_year,
        d.toxic = true
    """
    for batch in tqdm(list(batches(drugs, BATCH_SIZE)), desc="loading drugs"):
        session.run(query_drugs, rows=batch.to_dict("records"))

    print(f"[targets] {len(targets)} liens médicament-cible (mécanisme officiel)")
    query_targets = """
    UNWIND $rows AS row
    MATCH (d:Drug {chembl_id: row.chembl_id})
    MERGE (p:Protein {name: row.target_gene_symbol})
    MERGE (d)-[r:TARGETS]->(p)
    SET r.action_type = row.action_type,
        r.source = 'mechanism'
    """
    for batch in tqdm(list(batches(targets, BATCH_SIZE)), desc="loading targets (mechanism)"):
        session.run(query_targets, rows=batch.to_dict("records"))

    offtarget_path = PROC_DIR / "drug_targets_offtarget.csv"
    if offtarget_path.exists():
        offtargets = pd.read_csv(offtarget_path)
        print(f"[offtargets] {len(offtargets)} liens médicament-cible (bioactivité off-target)")
        query_offtargets = """
        UNWIND $rows AS row
        MATCH (d:Drug {chembl_id: row.chembl_id})
        MERGE (p:Protein {name: row.target_gene_symbol})
        MERGE (d)-[r:TARGETS]->(p)
        SET r.source = coalesce(r.source, 'bioactivity'),
            r.pchembl_value = row.pchembl_value,
            r.standard_type = row.standard_type
        """
        for batch in tqdm(list(batches(offtargets, BATCH_SIZE)), desc="loading targets (off-target)"):
            session.run(query_offtargets, rows=batch.to_dict("records"))
    else:
        print("[offtargets] fichier absent, lance 02b_download_chembl_offtargets.py si tu veux les off-targets")


def create_constraints(session):
    session.run("CREATE CONSTRAINT protein_name IF NOT EXISTS FOR (p:Protein) REQUIRE p.name IS UNIQUE")
    session.run("CREATE CONSTRAINT drug_id IF NOT EXISTS FOR (d:Drug) REQUIRE d.chembl_id IS UNIQUE")


def load_mutation_frequency(session):
    freq_path = PROC_DIR / "gene_mutation_frequency.csv"
    if not freq_path.exists():
        print("[mutations] fichier absent, lance scripts/08_download_cbioportal_mutations.py si tu veux ces données")
        return
    df = pd.read_csv(freq_path)
    print(f"[mutations] {len(df)} fréquences de mutation à charger")
    query = """
    UNWIND $rows AS row
    MATCH (p:Protein {name: row.gene})
    SET p.mutation_frequency = row.mutation_frequency,
        p.n_mutated_samples = row.n_mutated_samples
    """
    for batch in tqdm(list(batches(df, BATCH_SIZE)), desc="loading mutation frequency"):
        session.run(query, rows=batch.to_dict("records"))


def main():
    driver = GraphDatabase.driver(URI, auth=(USER, PASSWORD))
    with driver.session() as session:
        create_constraints(session)
        load_ppi(session)
        load_drugs_and_targets(session)
        load_mutation_frequency(session)
    driver.close()
    print("[done] import terminé")


if __name__ == "__main__":
    main()
