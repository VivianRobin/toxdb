"""
Compare les réseaux "influencés" (RWR/PageRank multi-seeds) de plusieurs
médicaments toxiques, pour tester l'hypothèse globale de la thèse :
des arrangements IPP communs se retrouvent entre principes actifs toxiques
partageant un même phénotype de toxicité.

Pour chaque médicament :
  1. Récupère ses cibles (TARGETS, mécanisme + off-target confondus)
  2. Lance un PageRank personnalisé (GDS) depuis ces cibles comme seeds
  3. Garde le top-N

Puis compare toutes les paires :
  - Similarité de Jaccard sur le top-N
  - Protéines communes / spécifiques à chaque médicament
  - Position (rang + score) d'une liste de gènes de référence optionnelle
    (utile pour vérifier si un module biologique connu -- ex: le complexe
    cardiaque du QT long -- est retrouvé chez chaque médicament)

Usage :
  python scripts/06_compare_drug_networks.py "ASTEMIZOLE" "TERFENADINE" \\
      --top-n 50 \\
      --reference-genes KCNQ1,KCNE1,KCNE2,SCN5A,CACNA1C,AKAP9,RYR2,KCNJ2

  # 3 médicaments ou plus : toutes les paires sont comparées
  python scripts/06_compare_drug_networks.py "ASTEMIZOLE" "TERFENADINE" "CISAPRIDE"

Sortie :
  - Résumé imprimé dans la console
  - CSV détaillé : data/processed/comparisons/<drugA>_vs_<drugB>.csv
"""
import argparse
import itertools
import os
from pathlib import Path

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


def ensure_graph_projected(session):
    exists = session.run(
        "CALL gds.graph.exists($name) YIELD exists RETURN exists", name=GRAPH_NAME
    ).single()["exists"]
    if not exists:
        print(f"[gds] projection '{GRAPH_NAME}' absente -> création")
        session.run(
            """
            CALL gds.graph.project(
              $name, 'Protein',
              {INTERACTS_WITH: {orientation: 'UNDIRECTED', properties: 'score'}}
            )
            """,
            name=GRAPH_NAME,
        )


def get_drug_ranking(session, drug_name: str, top_n: int) -> pd.DataFrame:
    """Retourne le classement PageRank complet (pas juste top_n) pour permettre
    de chercher le rang de gènes de référence même hors du top affiché."""
    query = """
    MATCH (d:Drug {name: $drug})-[:TARGETS]->(seed:Protein)
    WITH collect(seed) AS seeds
    CALL gds.pageRank.stream($graph, {
      sourceNodes: seeds,
      relationshipWeightProperty: 'score',
      dampingFactor: 0.85
    })
    YIELD nodeId, score
    RETURN gds.util.asNode(nodeId).name AS protein, score
    ORDER BY score DESC
    """
    records = session.run(query, drug=drug_name, graph=GRAPH_NAME).data()
    if not records:
        raise ValueError(
            f"Aucun résultat pour '{drug_name}'. Vérifie que le médicament existe "
            f"(MATCH (d:Drug) WHERE d.name CONTAINS '...' RETURN d.name) et qu'il a des cibles."
        )
    df = pd.DataFrame(records)
    df["rank"] = range(1, len(df) + 1)
    return df


def compare_pair(name_a: str, df_a: pd.DataFrame, name_b: str, df_b: pd.DataFrame, top_n: int, reference_genes: list[str]):
    top_a = set(df_a.head(top_n)["protein"])
    top_b = set(df_b.head(top_n)["protein"])

    intersection = top_a & top_b
    union = top_a | top_b
    jaccard = len(intersection) / len(union) if union else 0.0

    print(f"\n{'=' * 60}")
    print(f"{name_a}  vs  {name_b}   (top {top_n})")
    print(f"{'=' * 60}")
    print(f"Jaccard : {jaccard:.3f}  ({len(intersection)} communes / {len(union)} au total)")
    print(f"Communes ({len(intersection)}): {sorted(intersection)}")
    print(f"Spécifiques à {name_a} ({len(top_a - top_b)}): {sorted(top_a - top_b)}")
    print(f"Spécifiques à {name_b} ({len(top_b - top_a)}): {sorted(top_b - top_a)}")

    if reference_genes:
        print(f"\nGènes de référence (module biologique testé) :")
        rank_a = df_a.set_index("protein")["rank"].to_dict()
        rank_b = df_b.set_index("protein")["rank"].to_dict()
        n_total = len(df_a)
        for gene in reference_genes:
            ra = rank_a.get(gene, None)
            rb = rank_b.get(gene, None)
            ra_str = f"rang {ra}/{n_total}" if ra else "absent"
            rb_str = f"rang {rb}/{n_total}" if rb else "absent"
            print(f"  {gene:10s} | {name_a}: {ra_str:20s} | {name_b}: {rb_str}")

    # Sauvegarde détaillée
    COMPARISON_DIR.mkdir(parents=True, exist_ok=True)
    out = pd.DataFrame(
        {
            "protein": sorted(union),
        }
    )
    out[f"rank_{name_a}"] = out["protein"].map(df_a.set_index("protein")["rank"])
    out[f"score_{name_a}"] = out["protein"].map(df_a.set_index("protein")["score"])
    out[f"rank_{name_b}"] = out["protein"].map(df_b.set_index("protein")["rank"])
    out[f"score_{name_b}"] = out["protein"].map(df_b.set_index("protein")["score"])
    out["in_both_top_n"] = out["protein"].isin(intersection)
    out = out.sort_values(f"rank_{name_a}")

    safe_a = name_a.replace(" ", "_")
    safe_b = name_b.replace(" ", "_")
    out_path = COMPARISON_DIR / f"{safe_a}_vs_{safe_b}.csv"
    out.to_csv(out_path, index=False)
    print(f"\n[saved] {out_path}")

    return jaccard


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("drugs", nargs="+", help="Noms des médicaments (>=2), tels que stockés dans Neo4j (d.name)")
    parser.add_argument("--top-n", type=int, default=50, help="Taille du top-N comparé (défaut: 50)")
    parser.add_argument(
        "--reference-genes",
        type=str,
        default="",
        help="Liste de gènes séparés par des virgules, ex: KCNQ1,KCNE1,SCN5A,CACNA1C,AKAP9,RYR2",
    )
    args = parser.parse_args()

    if len(args.drugs) < 2:
        parser.error("Il faut au moins 2 médicaments à comparer")

    reference_genes = [g.strip() for g in args.reference_genes.split(",") if g.strip()]

    driver = GraphDatabase.driver(URI, auth=(USER, PASSWORD))
    with driver.session() as session:
        ensure_graph_projected(session)

        rankings = {}
        for drug in args.drugs:
            print(f"[pagerank] {drug}...")
            rankings[drug] = get_drug_ranking(session, drug, args.top_n)

    driver.close()

    jaccards = []
    for name_a, name_b in itertools.combinations(args.drugs, 2):
        j = compare_pair(name_a, rankings[name_a], name_b, rankings[name_b], args.top_n, reference_genes)
        jaccards.append({"drug_a": name_a, "drug_b": name_b, "jaccard": j})

    if len(args.drugs) > 2:
        print(f"\n{'=' * 60}\nRésumé toutes paires\n{'=' * 60}")
        summary = pd.DataFrame(jaccards).sort_values("jaccard", ascending=False)
        print(summary.to_string(index=False))
        summary.to_csv(COMPARISON_DIR / "summary_all_pairs.csv", index=False)


if __name__ == "__main__":
    main()
