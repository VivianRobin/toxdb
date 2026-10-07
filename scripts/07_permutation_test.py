"""
Test de permutation pour évaluer la significativité statistique du Jaccard
observé entre les réseaux PageRank de deux médicaments toxiques.

Modèle nul : à chaque permutation, on remplace les vraies cibles (seeds) de
chaque médicament par un ensemble aléatoire de MÊME TAILLE et de DEGRÉ DE
CONNECTIVITÉ SIMILAIRE (même bin de degré dans le réseau IPP), puis on
relance le PageRank et on recalcule le Jaccard du top-N.

Pourquoi apparier sur le degré : les cibles pharmacologiques réelles sont
souvent des protéines "hub" très connectées (récepteurs de signalisation).
Un tirage aléatoire uniforme donnerait des seeds peu connectées, un Jaccard
nul artificiellement bas, et donc une significativité gonflée de façon
trompeuse. Le modèle nul doit contrôler pour ce biais structurel.

p-value = fraction des permutations où le Jaccard nul >= Jaccard observé
          (test unilatéral : on teste si l'overlap réel est ANORMALEMENT ÉLEVÉ)

Usage :
  python scripts/07_permutation_test.py "ASTEMIZOLE" "TERFENADINE" \\
      --top-n 50 --n-permutations 200 --seed 42

Sortie :
  - p-value + z-score affichés en console
  - CSV des valeurs nulles : data/processed/comparisons/permtest_<A>_vs_<B>.csv
  - Histogramme (si matplotlib dispo) : data/processed/comparisons/permtest_<A>_vs_<B>.png
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


def get_all_degrees(session) -> pd.DataFrame:
    """Degré de chaque protéine (nb d'interactions IPP), pour l'appariement."""
    query = """
    MATCH (p:Protein)
    OPTIONAL MATCH (p)-[r:INTERACTS_WITH]-()
    RETURN p.name AS protein, count(r) AS degree
    """
    df = pd.DataFrame(session.run(query).data())
    return df


def build_degree_bins(degree_df: pd.DataFrame, n_bins: int) -> dict:
    """Retourne {protein_name: bin_id} et {bin_id: [protein_names]}."""
    # qcut sur le rang plutôt que la valeur brute -> bins de taille égale
    # même si la distribution de degré est très asymétrique (typique en PPI)
    degree_df = degree_df.copy()
    degree_df["bin"] = pd.qcut(degree_df["degree"].rank(method="first"), q=n_bins, labels=False)

    bin_to_proteins = defaultdict(list)
    protein_to_bin = {}
    for _, row in degree_df.iterrows():
        bin_to_proteins[row["bin"]].append(row["protein"])
        protein_to_bin[row["protein"]] = row["bin"]
    return protein_to_bin, bin_to_proteins


def get_real_seeds(session, drug_name: str) -> list[str]:
    query = "MATCH (:Drug {name: $drug})-[:TARGETS]->(p:Protein) RETURN DISTINCT p.name AS protein"
    result = session.run(query, drug=drug_name).data()
    if not result:
        raise ValueError(f"Aucune cible trouvée pour '{drug_name}'")
    return [r["protein"] for r in result]


def pagerank_top_n(session, seed_names: list[str], top_n: int) -> set:
    query = """
    MATCH (p:Protein) WHERE p.name IN $seed_names
    WITH collect(p) AS seeds
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
    result = session.run(query, seed_names=seed_names, graph=GRAPH_NAME, top_n=top_n).data()
    return {r["protein"] for r in result}


def sample_degree_matched(real_seeds: list[str], protein_to_bin: dict, bin_to_proteins: dict, rng: random.Random) -> list[str]:
    sampled = []
    for seed in real_seeds:
        b = protein_to_bin.get(seed)
        candidates = bin_to_proteins.get(b, [])
        if not candidates:
            continue
        sampled.append(rng.choice(candidates))
    return sampled


def jaccard(a: set, b: set) -> float:
    union = a | b
    if not union:
        return 0.0
    return len(a & b) / len(union)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("drug_a")
    parser.add_argument("drug_b")
    parser.add_argument("--top-n", type=int, default=50)
    parser.add_argument("--n-permutations", type=int, default=200)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)

    driver = GraphDatabase.driver(URI, auth=(USER, PASSWORD))
    with driver.session() as session:
        ensure_graph_projected(session)

        print("[setup] calcul des degrés pour l'appariement...")
        degree_df = get_all_degrees(session)
        protein_to_bin, bin_to_proteins = build_degree_bins(degree_df, N_DEGREE_BINS)

        print(f"[observed] {args.drug_a} vs {args.drug_b}")
        seeds_a = get_real_seeds(session, args.drug_a)
        seeds_b = get_real_seeds(session, args.drug_b)
        print(f"  {args.drug_a}: {len(seeds_a)} cibles réelles")
        print(f"  {args.drug_b}: {len(seeds_b)} cibles réelles")

        top_a_obs = pagerank_top_n(session, seeds_a, args.top_n)
        top_b_obs = pagerank_top_n(session, seeds_b, args.top_n)
        observed_jaccard = jaccard(top_a_obs, top_b_obs)
        print(f"  Jaccard observé : {observed_jaccard:.4f}")

        print(f"\n[permutations] {args.n_permutations} tirages, seeds appariées par degré...")
        null_jaccards = []
        for i in range(args.n_permutations):
            rand_seeds_a = sample_degree_matched(seeds_a, protein_to_bin, bin_to_proteins, rng)
            rand_seeds_b = sample_degree_matched(seeds_b, protein_to_bin, bin_to_proteins, rng)

            top_a_null = pagerank_top_n(session, rand_seeds_a, args.top_n)
            top_b_null = pagerank_top_n(session, rand_seeds_b, args.top_n)
            null_jaccards.append(jaccard(top_a_null, top_b_null))

            if (i + 1) % 50 == 0:
                print(f"  {i + 1}/{args.n_permutations}...")

    driver.close()

    null_jaccards = np.array(null_jaccards)
    p_value = (np.sum(null_jaccards >= observed_jaccard) + 1) / (args.n_permutations + 1)
    z_score = (observed_jaccard - null_jaccards.mean()) / null_jaccards.std() if null_jaccards.std() > 0 else float("inf")

    print(f"\n{'=' * 60}")
    print("RÉSULTAT")
    print(f"{'=' * 60}")
    print(f"Jaccard observé          : {observed_jaccard:.4f}")
    print(f"Jaccard nul (moyenne)    : {null_jaccards.mean():.4f}  (écart-type: {null_jaccards.std():.4f})")
    print(f"Jaccard nul (min-max)    : {null_jaccards.min():.4f} - {null_jaccards.max():.4f}")
    print(f"z-score                  : {z_score:.2f}")
    print(f"p-value                  : {p_value:.4f}  ({'significatif' if p_value < 0.05 else 'NON significatif'} à alpha=0.05)")

    COMPARISON_DIR.mkdir(parents=True, exist_ok=True)
    safe_a = args.drug_a.replace(" ", "_")
    safe_b = args.drug_b.replace(" ", "_")

    out_csv = COMPARISON_DIR / f"permtest_{safe_a}_vs_{safe_b}.csv"
    pd.DataFrame({"null_jaccard": null_jaccards}).to_csv(out_csv, index=False)
    print(f"\n[saved] {out_csv}")

    try:
        import matplotlib
        matplotlib.use("Agg")  # backend non-interactif, pas besoin d'affichage graphique
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(7, 4))
        ax.hist(null_jaccards, bins=30, color="#94a3b8", edgecolor="white", label="Distribution nulle\n(seeds aléatoires, appariées par degré)")
        ax.axvline(observed_jaccard, color="#dc2626", linewidth=2, label=f"Jaccard observé = {observed_jaccard:.3f}")
        ax.set_xlabel("Jaccard (top-{})".format(args.top_n))
        ax.set_ylabel("Nombre de permutations")
        ax.set_title(f"{args.drug_a} vs {args.drug_b}\np = {p_value:.4f}, z = {z_score:.2f}")
        ax.legend(fontsize=8)
        fig.tight_layout()
        out_png = COMPARISON_DIR / f"permtest_{safe_a}_vs_{safe_b}.png"
        fig.savefig(out_png, dpi=150)
        print(f"[saved] {out_png}")
    except ImportError:
        print("[info] matplotlib non installé, histogramme non généré (pip install matplotlib)")


if __name__ == "__main__":
    main()
