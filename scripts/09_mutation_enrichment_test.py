"""
Teste si un ensemble de gènes (l'arrangement IPP toxique commun, ou le
top-N PageRank d'un médicament) est enrichi en mutations pan-cancer (TCGA
PanCancer Atlas, via scripts/08_download_cbioportal_mutations.py) par
rapport au hasard.

Même logique que 07_permutation_test.py : modèle nul avec gènes aléatoires
appariés par degré de connectivité (pour ne pas confondre "hub proteins
sont plus mutées" avec un vrai enrichissement toxicité-mutation).

Hypothèse testée : les gènes de l'arrangement toxique ont une fréquence de
mutation pan-cancer MOYENNE plus élevée que des gènes aléatoires de même
degré dans le réseau IPP.

Usage :
  # Sur le top-N PageRank d'un médicament
  python scripts/09_mutation_enrichment_test.py --drug ASTEMIZOLE --top-n 50 --n-permutations 500

  # Sur une liste explicite (ex: l'arrangement cardiaque qu'on a identifié)
  python scripts/09_mutation_enrichment_test.py \\
      --genes KCNQ1,KCNE1,KCNE2,SCN5A,CACNA1C,AKAP9,RYR2,KCNJ2 \\
      --n-permutations 500

Prérequis : avoir chargé gene_mutation_frequency.csv dans Neo4j
            (scripts/03_load_neo4j.py après 08_download_cbioportal_mutations.py)
"""
import argparse
import os
import random
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from dotenv import load_dotenv
from neo4j import GraphDatabase

load_dotenv()

URI = os.getenv("NEO4J_URI", "bolt://localhost:7687")
USER = os.getenv("NEO4J_USER", "neo4j")
PASSWORD = os.getenv("NEO4J_PASSWORD", "toxkg-dev-pw")

PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"
COMPARISON_DIR = PROC_DIR / "comparisons"

GRAPH_NAME = "ppi"
N_DEGREE_BINS = 20


def ensure_graph_projected(session):
    exists = session.run(
        "CALL gds.graph.exists($name) YIELD exists RETURN exists", name=GRAPH_NAME
    ).single()["exists"]
    if not exists:
        session.run(
            """
            CALL gds.graph.project(
              $name, 'Protein',
              {INTERACTS_WITH: {orientation: 'UNDIRECTED', properties: 'score'}}
            )
            """,
            name=GRAPH_NAME,
        )


def get_all_degrees_and_mutfreq(session) -> pd.DataFrame:
    query = """
    MATCH (p:Protein)
    OPTIONAL MATCH (p)-[r:INTERACTS_WITH]-()
    WITH p, count(r) AS degree
    RETURN p.name AS protein, degree, p.mutation_frequency AS mutation_frequency
    """
    df = pd.DataFrame(session.run(query).data())
    n_missing = df["mutation_frequency"].isna().sum()
    if n_missing > 0:
        print(
            f"[warning] {n_missing}/{len(df)} protéines sans fréquence de mutation "
            f"(pas dans gene_mutation_frequency.csv) -- elles seront exclues du tirage aléatoire"
        )
    return df.dropna(subset=["mutation_frequency"])


def build_degree_bins(df: pd.DataFrame, n_bins: int) -> tuple[dict, dict]:
    df = df.copy()
    df["bin"] = pd.qcut(df["degree"].rank(method="first"), q=n_bins, labels=False)
    bin_to_proteins = defaultdict(list)
    protein_to_bin = {}
    for _, row in df.iterrows():
        bin_to_proteins[row["bin"]].append(row["protein"])
        protein_to_bin[row["protein"]] = row["bin"]
    return protein_to_bin, bin_to_proteins


def get_drug_top_n(session, drug_name: str, top_n: int) -> list[str]:
    query = """
    MATCH (d:Drug {name: $drug})-[:TARGETS]->(seed:Protein)
    WITH collect(seed) AS seeds
    CALL gds.pageRank.stream($graph, {
      sourceNodes: seeds,
      relationshipWeightProperty: 'score',
      dampingFactor: 0.85
    })
    YIELD nodeId, score
    RETURN gds.util.asNode(nodeId).name AS protein
    ORDER BY score DESC
    LIMIT $top_n
    """
    result = session.run(query, drug=drug_name, graph=GRAPH_NAME, top_n=top_n).data()
    return [r["protein"] for r in result]


def sample_degree_matched(genes: list[str], protein_to_bin: dict, bin_to_proteins: dict, rng: random.Random) -> list[str]:
    sampled = []
    for g in genes:
        b = protein_to_bin.get(g)
        candidates = bin_to_proteins.get(b, [])
        if candidates:
            sampled.append(rng.choice(candidates))
    return sampled


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--drug", type=str, help="Nom du médicament (utilise son top-N PageRank)")
    src.add_argument("--genes", type=str, help="Liste de gènes séparés par virgules")
    parser.add_argument("--top-n", type=int, default=50, help="Si --drug: taille du top-N (défaut: 50)")
    parser.add_argument("--n-permutations", type=int, default=500)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)

    driver = GraphDatabase.driver(URI, auth=(USER, PASSWORD))
    with driver.session() as session:
        ensure_graph_projected(session)

        print("[setup] degrés + fréquences de mutation...")
        full_df = get_all_degrees_and_mutfreq(session)
        protein_to_bin, bin_to_proteins = build_degree_bins(full_df, N_DEGREE_BINS)
        mutfreq_lookup = dict(zip(full_df["protein"], full_df["mutation_frequency"]))

        if args.drug:
            print(f"[genes] top-{args.top_n} PageRank de {args.drug}")
            target_genes = get_drug_top_n(session, args.drug, args.top_n)
            label = f"{args.drug}_top{args.top_n}"
        else:
            target_genes = [g.strip() for g in args.genes.split(",") if g.strip()]
            label = "custom_geneset"

    driver.close()

    target_genes_with_freq = [g for g in target_genes if g in mutfreq_lookup]
    n_missing = len(target_genes) - len(target_genes_with_freq)
    if n_missing:
        print(f"[warning] {n_missing} gène(s) de l'ensemble testé sans fréquence de mutation, exclus")

    if not target_genes_with_freq:
        raise SystemExit("Aucun gène testable (fréquence de mutation manquante partout). Lance script 08 d'abord.")

    observed_freqs = [mutfreq_lookup[g] for g in target_genes_with_freq]
    observed_mean = float(np.mean(observed_freqs))

    print(f"\n{len(target_genes_with_freq)} gènes testés")
    print(f"Fréquence de mutation moyenne observée : {observed_mean:.4f}")

    print(f"\n[permutations] {args.n_permutations} tirages, gènes appariés par degré...")
    null_means = []
    for i in range(args.n_permutations):
        rand_genes = sample_degree_matched(target_genes_with_freq, protein_to_bin, bin_to_proteins, rng)
        rand_freqs = [mutfreq_lookup[g] for g in rand_genes if g in mutfreq_lookup]
        if rand_freqs:
            null_means.append(np.mean(rand_freqs))
        if (i + 1) % 100 == 0:
            print(f"  {i + 1}/{args.n_permutations}...")

    null_means = np.array(null_means)
    p_value = (np.sum(null_means >= observed_mean) + 1) / (len(null_means) + 1)
    z_score = (observed_mean - null_means.mean()) / null_means.std() if null_means.std() > 0 else float("inf")

    print(f"\n{'=' * 60}")
    print(f"RÉSULTAT -- {label}")
    print(f"{'=' * 60}")
    print(f"Fréquence de mutation observée (moyenne) : {observed_mean:.4f}")
    print(f"Fréquence nulle (moyenne des permutations): {null_means.mean():.4f}  (écart-type: {null_means.std():.4f})")
    print(f"z-score                                   : {z_score:.2f}")
    print(f"p-value                                   : {p_value:.4f}  ({'significatif' if p_value < 0.05 else 'NON significatif'} à alpha=0.05)")

    print("\nDétail par gène :")
    detail = pd.DataFrame(
        {"gene": target_genes_with_freq, "mutation_frequency": observed_freqs}
    ).sort_values("mutation_frequency", ascending=False)
    print(detail.to_string(index=False))

    COMPARISON_DIR.mkdir(parents=True, exist_ok=True)
    out_path = COMPARISON_DIR / f"mutation_enrichment_{label}.csv"
    pd.DataFrame({"null_mean_mutation_freq": null_means}).to_csv(out_path, index=False)
    detail.to_csv(COMPARISON_DIR / f"mutation_enrichment_{label}_detail.csv", index=False)
    print(f"\n[saved] {out_path}")


if __name__ == "__main__":
    main()
