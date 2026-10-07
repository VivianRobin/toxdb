from pathlib import Path
import pandas as pd
import igraph as ig

PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"
COMPARISON_DIR = PROC_DIR / "comparisons"
MODULE_DIR = PROC_DIR / "modules"
MODULE_DIR.mkdir(exist_ok=True, parents=True)

def main():
    print("--- 1. Chargement des données ---")
    ppi_edges = pd.read_csv(PROC_DIR / "ppi_edges.csv")
    crossref = pd.read_csv(COMPARISON_DIR / "motif_mutation_crossref.csv")
    
    # Filtrage optionnel (ex: score >= 400)
    if "score" in ppi_edges.columns:
        initial_edges = len(ppi_edges)
        ppi_edges = ppi_edges[ppi_edges["score"] >= 400]
        print(f"Filtrage des arêtes (score >= 400) : {len(ppi_edges)} / {initial_edges} arêtes conservées.")

    print("\n--- 2. Construction du graphe avec igraph (C-optimized) ---")
    # Correction : Suppression de vertex_name_attr (géré automatiquement par igraph)
    G = ig.Graph.DataFrame(ppi_edges, directed=False,use_vids=False)
    print(f"Réseau global chargé : {G.vcount()} nœuds, {G.ecount()} arêtes.")

    print("\n--- 3. Détection de communautés (Algorithme de Louvain / Multilevel) ---")
    weights = "score" if "score" in ppi_edges.columns else None
    communities = G.community_multilevel(weights=weights)
    
    print(f"{len(communities)} modules/communautés détectés instantanément.")

    priority_map = crossref.set_index("gene")["priority_score"].to_dict()
    class_map = crossref.set_index("gene")["warning_class"].to_dict()

    module_summary = []

    for idx, comm_indices in enumerate(communities):
        comm_nodes = [G.vs[node]["name"] for node in comm_indices]
        
        if len(comm_nodes) < 5:
            continue
        
        subG = G.induced_subgraph(comm_indices)
        density = subG.density()
        
        comm_tox_genes = [g for g in comm_nodes if g in priority_map]
        
        if comm_tox_genes:
            classes_in_mod = [class_map[g] for g in comm_tox_genes]
            dominant_class = max(set(classes_in_mod), key=classes_in_mod.count)
            avg_priority = sum(priority_map[g] for g in comm_tox_genes) / len(comm_tox_genes)
            
            module_summary.append({
                "module_id": idx,
                "size": len(comm_nodes),
                "density": round(density, 3),
                "dominant_toxicity_class": dominant_class,
                "num_tox_genes": len(comm_tox_genes),
                "avg_priority_score": round(avg_priority, 3),
                "genes": ", ".join(comm_nodes[:10])
            })

    df_modules = pd.DataFrame(module_summary)
    if not df_modules.empty:
        df_modules = df_modules.sort_values(by="avg_priority_score", ascending=False)
        out_path = MODULE_DIR / "dense_modules_summary.csv"
        df_modules.to_csv(out_path, index=False)
        print(f"\n[saved] Résumé des modules enregistré dans {out_path}")
        print(df_modules[["module_id", "size", "density", "dominant_toxicity_class", "num_tox_genes", "avg_priority_score"]].head(10).to_string(index=False))
    else:
        print("[info] Aucun module dense significatif croisé avec les gènes toxiques.")

if __name__ == "__main__":
    main()
