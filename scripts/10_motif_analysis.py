"""
Recherche de motifs de réseau (cliques maximales dans le voisinage PageRank
des cibles) pour L'ENSEMBLE des médicaments toxiques, et test si les
médicaments qui partagent un même phénotype de toxicité (d.warning_class,
ex: cardiotoxicity, hepatotoxicity...) partagent significativement plus de
"motif nodes" que des regroupements aléatoires de médicaments.

Pourquoi cette approche (et pas un simple 07_permutation_test.py répété
sur toutes les paires) : avec ~170 médicaments ayant des cibles connues,
il y a ~14000 paires possibles. Relancer un PageRank + un modèle nul par
paire (comme 07) est infaisable à cette échelle. On sépare donc le calcul
en deux temps :

  1. UNE SEULE fois par médicament (pas par paire) :
       a. PageRank personnalisé (GDS) depuis ses cibles (TARGETS, mécanisme
          + off-target confondus) comme seeds
       b. Sous-graphe IPP induit sur le top-N du classement PageRank
       c. Cliques maximales (networkx.find_cliques) dans ce sous-graphe
       d. "motif nodes" = union des protéines appartenant à >= 1 clique
          de taille >= --min-clique-size (les cliques isolées de taille 2
          sont juste des arêtes -- pas d'intérêt structurel)

  2. Similarité de Jaccard (motif nodes) pour TOUTES les paires de
     médicaments -- calcul en mémoire pur (numpy), aucune requête Neo4j
     supplémentaire, donc instantané même sur ~14000 paires.

  3. Test de significativité par PERMUTATION D'ÉTIQUETTES : on mélange
     aléatoirement les warning_class entre médicaments (la topologie
     réelle de chaque médicament est conservée, seule l'affectation
     phénotype <-> médicament est randomisée), et on compare le Jaccard
     moyen intra-classe observé à la distribution nulle. C'est un modèle
     nul différent de celui de 07/09 (qui rééchantillonne des seeds
     appariées par degré) : ici on teste "est-ce que le REGROUPEMENT par
     phénotype capture plus de similarité structurelle que n'importe quel
     regroupement aléatoire de même taille parmi les médicaments réels ?"

Usage :
  python scripts/10_motif_analysis.py
  python scripts/10_motif_analysis.py --top-n 40 --min-clique-size 3 \\
      --min-targets 1 --min-drugs-per-class 3 --n-permutations 1000

Sorties (data/processed/comparisons/) :
  - drug_motif_nodes.csv           1 ligne / médicament (cliques, motif nodes)
  - drug_pair_motif_jaccard.csv    toutes les paires (jaccard, same_class)
  - motif_enrichment_by_class.csv  résumé par warning_class (obs/null/p/z)
"""
import argparse
import itertools
import os
import random
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from dotenv import load_dotenv
from neo4j import GraphDatabase
from tqdm import tqdm

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


def get_drugs_with_targets(session, min_targets: int) -> pd.DataFrame:
    """Tous les médicaments toxiques ayant >= min_targets cibles connues
    (mécanisme officiel + off-target bioactivité confondus)."""
    query = """
    MATCH (d:Drug)-[:TARGETS]->(p:Protein)
    WITH d, count(DISTINCT p) AS n_targets
    WHERE n_targets >= $min_targets
    RETURN d.chembl_id AS chembl_id, d.name AS name, d.warning_class AS warning_class,
           n_targets
    ORDER BY d.name
    """
    df = pd.DataFrame(session.run(query, min_targets=min_targets).data())
    return df


def get_seeds(session, chembl_id: str) -> list[str]:
    query = "MATCH (:Drug {chembl_id: $id})-[:TARGETS]->(p:Protein) RETURN DISTINCT p.name AS protein"
    return [r["protein"] for r in session.run(query, id=chembl_id).data()]


def pagerank_top_n(session, seed_names: list[str], top_n: int) -> list[str]:
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
    return [r["protein"] for r in result]


def induced_subgraph_edges(session, node_names: list[str]) -> nx.Graph:
    """Arêtes IPP RÉELLES entre les nœuds du top-N (nécessaire pour que
    les cliques aient un sens -- le top-N PageRank n'est qu'une liste,
    pas un graphe)."""
    query = """
    MATCH (a:Protein)-[r:INTERACTS_WITH]-(b:Protein)
    WHERE a.name IN $nodes AND b.name IN $nodes AND a.name < b.name
    RETURN a.name AS source, b.name AS target, r.score AS score
    """
    G = nx.Graph()
    G.add_nodes_from(node_names)
    for record in session.run(query, nodes=node_names).data():
        G.add_edge(record["source"], record["target"], score=record["score"])
    return G


def motif_nodes_for_drug(session, chembl_id: str, top_n: int, min_clique_size: int):
    seeds = get_seeds(session, chembl_id)
    top_nodes = pagerank_top_n(session, seeds, top_n)
    if len(top_nodes) < min_clique_size:
        return set(), 0, 0

    G = induced_subgraph_edges(session, top_nodes)
    cliques = [c for c in nx.find_cliques(G) if len(c) >= min_clique_size]
    motif_nodes = set().union(*cliques) if cliques else set()
    max_clique_size = max((len(c) for c in cliques), default=0)
    return motif_nodes, len(cliques), max_clique_size


def jaccard(a: set, b: set) -> float:
    union = a | b
    if not union:
        return 0.0
    return len(a & b) / len(union)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--top-n", type=int, default=40, help="Taille du voisinage PageRank par médicament (défaut: 40)")
    parser.add_argument("--min-clique-size", type=int, default=3, help="Taille minimale d'une clique pour compter comme motif (défaut: 3)")
    parser.add_argument("--min-targets", type=int, default=1, help="Nb minimal de cibles connues pour inclure un médicament (défaut: 1)")
    parser.add_argument("--min-drugs-per-class", type=int, default=3, help="Nb minimal de médicaments par warning_class pour tester cette classe (défaut: 3)")
    parser.add_argument("--n-permutations", type=int, default=1000, help="Nb de permutations d'étiquettes (défaut: 1000)")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    COMPARISON_DIR.mkdir(parents=True, exist_ok=True)

    driver = GraphDatabase.driver(URI, auth=(USER, PASSWORD))
    with driver.session() as session:
        ensure_graph_projected(session)

        print(f"[setup] médicaments avec >= {args.min_targets} cible(s)...")
        drugs = get_drugs_with_targets(session, args.min_targets)
        print(f"  -> {len(drugs)} médicaments")

        print(f"\n[1/3] PageRank + cliques par médicament (top-{args.top_n}, clique >= {args.min_clique_size})...")
        motif_nodes_list = []
        n_cliques_list = []
        max_clique_list = []
        for chembl_id in tqdm(drugs["chembl_id"], desc="motifs"):
            nodes, n_cliques, max_size = motif_nodes_for_drug(session, chembl_id, args.top_n, args.min_clique_size)
            motif_nodes_list.append(nodes)
            n_cliques_list.append(n_cliques)
            max_clique_list.append(max_size)

    driver.close()

    drugs = drugs.copy()
    drugs["n_cliques"] = n_cliques_list
    drugs["max_clique_size"] = max_clique_list
    drugs["n_motif_nodes"] = [len(s) for s in motif_nodes_list]

    out_drugs = drugs.copy()
    out_drugs["motif_nodes"] = ["|".join(sorted(s)) for s in motif_nodes_list]
    out_path = COMPARISON_DIR / "drug_motif_nodes.csv"
    out_drugs.to_csv(out_path, index=False)
    print(f"[saved] {out_path}")
    print(
        f"  médiane cliques/médicament: {drugs['n_cliques'].median():.0f} | "
        f"médiane motif nodes/médicament: {drugs['n_motif_nodes'].median():.0f}"
    )

    # --- 2. Jaccard toutes paires (en mémoire) ---
    print("\n[2/3] Jaccard (motif nodes) sur toutes les paires...")
    n = len(drugs)
    idx_a, idx_b = np.triu_indices(n, k=1)
    jac = np.array([jaccard(motif_nodes_list[i], motif_nodes_list[j]) for i, j in zip(idx_a, idx_b)])

    names = drugs["name"].to_numpy()
    classes = drugs["warning_class"].to_numpy()

    pairs_df = pd.DataFrame(
        {
            "drug_a": names[idx_a],
            "drug_b": names[idx_b],
            "class_a": classes[idx_a],
            "class_b": classes[idx_b],
            "jaccard_motif_nodes": jac,
        }
    )
    pairs_df["same_class"] = (pairs_df["class_a"] == pairs_df["class_b"]) & pairs_df["class_a"].notna()
    pairs_path = COMPARISON_DIR / "drug_pair_motif_jaccard.csv"
    pairs_df.sort_values("jaccard_motif_nodes", ascending=False).to_csv(pairs_path, index=False)
    print(f"[saved] {pairs_path}  ({len(pairs_df)} paires)")

    # --- 3. Test par permutation d'étiquettes, restreint aux classes assez peuplées ---
    print(f"\n[3/3] Test d'enrichissement par phénotype (permutation d'étiquettes, {args.n_permutations} tirages)...")

    class_counts = pd.Series(classes).value_counts(dropna=True)
    eligible_classes = class_counts[class_counts >= args.min_drugs_per_class].index.tolist()
    eligible_mask = np.isin(classes, eligible_classes)

    if not eligible_classes:
        print(f"  Aucune classe avec >= {args.min_drugs_per_class} médicaments -- rien à tester.")
        return

    # On restreint la matrice de paires aux seuls médicaments dont la classe
    # est testable, pour que le fond de permutation soit homogène.
    sub_idx = np.where(eligible_mask)[0]
    sub_labels = classes[sub_idx]
    sub_i, sub_j = np.triu_indices(len(sub_idx), k=1)

    # matrice complète -> sous-ensemble : reconstruit jac pour les paires sub
    full_jac_matrix = np.zeros((n, n))
    full_jac_matrix[idx_a, idx_b] = jac
    full_jac_matrix[idx_b, idx_a] = jac
    sub_jac = full_jac_matrix[np.ix_(sub_idx, sub_idx)]
    pair_jac = sub_jac[sub_i, sub_j]

    def mean_intra_class(labels: np.ndarray) -> dict:
        la, lb = labels[sub_i], labels[sub_j]
        same = la == lb
        out = {}
        for c in eligible_classes:
            m = same & (la == c)
            out[c] = pair_jac[m].mean() if m.any() else np.nan
        return out

    observed = mean_intra_class(sub_labels)

    null_dist = {c: np.empty(args.n_permutations) for c in eligible_classes}
    perm_labels = sub_labels.copy()
    for p in range(args.n_permutations):
        rng.shuffle(perm_labels)
        res = mean_intra_class(perm_labels)
        for c in eligible_classes:
            null_dist[c][p] = res[c]

    rows = []
    for c in eligible_classes:
        obs = observed[c]
        null = null_dist[c]
        null_valid = null[~np.isnan(null)]
        p_value = (np.sum(null_valid >= obs) + 1) / (len(null_valid) + 1) if len(null_valid) else np.nan
        z = (obs - null_valid.mean()) / null_valid.std() if len(null_valid) and null_valid.std() > 0 else np.nan
        rows.append(
            {
                "warning_class": c,
                "n_drugs": int(class_counts[c]),
                "jaccard_intra_class_observed": obs,
                "jaccard_intra_class_null_mean": null_valid.mean() if len(null_valid) else np.nan,
                "jaccard_intra_class_null_std": null_valid.std() if len(null_valid) else np.nan,
                "z_score": z,
                "p_value": p_value,
                "significant_0.05": bool(p_value < 0.05) if not np.isnan(p_value) else False,
            }
        )

    summary = pd.DataFrame(rows).sort_values("z_score", ascending=False)
    summary_path = COMPARISON_DIR / "motif_enrichment_by_class.csv"
    summary.to_csv(summary_path, index=False)

    print(f"\n{'=' * 78}")
    print("RÉSULTAT -- enrichissement en motifs de réseau par phénotype de toxicité")
    print(f"{'=' * 78}")
    print(summary.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print(f"\n[saved] {summary_path}")


if __name__ == "__main__":
    main()
