from pathlib import Path
import pandas as pd
import igraph as ig
import numpy as np

PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"
COMPARISON_DIR = PROC_DIR / "comparisons"
RWR_DIR = PROC_DIR / "rwr_results"
RWR_DIR.mkdir(exist_ok=True, parents=True)

def main():
    print("--- 1. Chargement des données pour le RWR ---")
    ppi_edges = pd.read_csv(PROC_DIR / "ppi_edges.csv")
    crossref = pd.read_csv(COMPARISON_DIR / "motif_mutation_crossref.csv")
    
    # Filtrer optionnellement par score
    if "score" in ppi_edges.columns:
        ppi_edges = ppi_edges[ppi_edges["score"] >= 400]

    print("\n--- 2. Construction du graphe igraph ---")
    G = ig.Graph.DataFrame(ppi_edges, directed=False, use_vids=False)
    print(f"Réseau chargé : {G.vcount()} nœuds, {G.ecount()} arêtes.")

    # Obtenir la liste de tous les noms de nœuds du graphe pour aligner les vecteurs
    all_nodes = [v["name"] for v in G.vs]
    node_to_idx = {name: i for i, name in enumerate(all_nodes)}

    # Dictionnaire des priorités
    priority_map = crossref.set_index("gene")["priority_score"].to_dict()
    classes = crossref["warning_class"].unique()

    weights = "score" if "score" in ppi_edges.columns else None

    print("\n--- 3. Lancement de la Marche Aléatoire avec Redémarrage (RWR) par classe ---")
    
    all_rwr_results = []

    for cls in classes:
        print(f"\nTraitement de la classe : {cls}")
        
        # Sélectionner les graines (seed genes) pour cette classe (ex: Top 20 gènes priorisés)
        df_cls = crossref[crossref["warning_class"] == cls].sort_values(by="priority_score", ascending=False)
        seed_genes = df_cls["gene"].head(20).tolist()
        
        # Filtrer pour s'assurer qu'ils sont bien dans le graphe igraph
        seed_genes = [g for g in seed_genes if g in node_to_idx]
        
        if not seed_genes:
            print(f"[skip] Pas assez de graines valides pour {cls}")
            continue

        # Construire le vecteur de redémarrage (reset vector)
        reset_vector = np.zeros(G.vcount())
        for g in seed_genes:
            reset_vector[node_to_idx[g]] = 1.0
        
        # Normaliser pour que la somme fasse 1 (requis par igraph)
        if reset_vector.sum() > 0:
            reset_vector /= reset_vector.sum()

        # Exécuter le PageRank personnalisé (qui équivaut mathématiquement au RWR)
        # damping=0.85 signifie 15% de chance de redémarrer sur les graines à chaque pas
        pr_scores = G.personalized_pagerank(
            reset=reset_vector.tolist(),
            damping=0.85,
            weights=weights
        )

        # Créer un DataFrame des scores pour tous les nœuds du réseau
        df_pr = pd.DataFrame({
            "gene": all_nodes,
            "rwr_score": pr_scores,
            "warning_class": cls
        })

        # Marquer s'il s'agit d'une graine initiale ou d'une découverte (gêne fantôme)
        df_pr["is_seed"] = df_pr["gene"].isin(seed_genes)

        # Trier par score RWR décroissant
        df_pr = df_pr.sort_values(by="rwr_score", ascending=False)
        
        # Afficher le Top 5 des gènes découverts (hors graines initiales)
        ghosts = df_pr[~df_pr["is_seed"]].head(5)
        print(f"Top 5 'gènes fantômes' découverts par propagation pour {cls} :")
        for _, row in ghosts.iterrows():
            print(f"  -> {row['gene']} (score RWR: {row['rwr_score']:.5f})")

        all_rwr_results.append(df_pr.head(50)) # Garder le top 50 par classe

    # Sauvegarder les résultats globaux
    final_rwr_df = pd.concat(all_rwr_results, ignore_index=True)
    out_path = RWR_DIR / "rwr_top_propagated_genes.csv"
    final_rwr_df.to_csv(out_path, index=False)
    print(f"\n[saved] Résultats de propagation RWR enregistrés dans {out_path}")

if __name__ == "__main__":
    main()
