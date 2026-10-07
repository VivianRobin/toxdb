import argparse
from pathlib import Path
import pandas as pd
import networkx as nx

# Force le backend non-interactif de matplotlib pour éviter l'erreur d'affichage (Tkinter/display)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"
COMPARISON_DIR = PROC_DIR / "comparisons"
VIZ_DIR = PROC_DIR / "visualizations"
VIZ_DIR.mkdir(exist_ok=True, parents=True)

def main():
    parser = argparse.ArgumentParser(description="Génère les graphes des sous-réseaux IPP par classe.")
    parser.add_argument("--top-n", type=int, default=15, help="Nb de gènes à inclure dans le réseau (defaut: 15)")
    args = parser.parse_args()

    # 1. Charger les données
    crossref = pd.read_csv(COMPARISON_DIR / "motif_mutation_crossref.csv")
    ppi_edges = pd.read_csv(PROC_DIR / "ppi_edges.csv")

    classes = crossref["warning_class"].unique()

    for cls in classes:
        # 2. Isoler les gènes du Top N pour cette classe
        df_cls = crossref[crossref["warning_class"] == cls].head(args.top_n)
        top_genes = set(df_cls["gene"])
        
        if not top_genes:
            continue

        # 3. Filtrer les interactions IPP : on ne garde que les arêtes entre nos top gènes
        sub_edges = ppi_edges[
            (ppi_edges["source"].isin(top_genes)) & 
            (ppi_edges["target"].isin(top_genes))
        ]

        # 4. Construire le graphe
        G = nx.from_pandas_edgelist(sub_edges, source="source", target="target")
        
        # Ajouter les gènes isolés (qui font partie du top N mais n'ont pas d'interaction directe)
        G.add_nodes_from(top_genes)

        # 5. Visualisation
        plt.figure(figsize=(10, 8))
        plt.title(f"Sous-réseau IPP : {cls.capitalize()} (Top {args.top_n})", fontsize=16)
        
        # Mapping des tailles de nœuds basé sur le priority_score
        scores = df_cls.set_index("gene")["priority_score"].to_dict()
        node_sizes = [scores.get(node, 0.5) * 1500 for node in G.nodes()]

        # Dessiner le réseau
        pos = nx.spring_layout(G, k=0.5, seed=42)
        nx.draw(
            G, pos, 
            with_labels=True, 
            node_size=node_sizes, 
            node_color="skyblue", 
            edge_color="gray", 
            font_size=10, 
            font_weight="bold",
            alpha=0.9
        )
        
        out_path = VIZ_DIR / f"network_{cls.replace(' ', '_')}.png"
        plt.savefig(out_path, format="png", dpi=300, bbox_inches="tight")
        plt.close()
        
        print(f"[saved] {out_path} ({G.number_of_nodes()} nœuds, {G.number_of_edges()} arêtes)")

if __name__ == "__main__":
    main()
