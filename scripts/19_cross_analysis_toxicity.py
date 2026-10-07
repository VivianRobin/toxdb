from pathlib import Path
import pandas as pd

PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"
CLINICAL_DIR = PROC_DIR / "cbioportal"
COMP_DIR = PROC_DIR / "comparisons"

def main():
    print("--- Croisement des altérations CNA avec les classes de toxicité ---")
    
    # 1. Chargement des données
    stats_path = CLINICAL_DIR / "top_altered_genes_summary.csv"
    crossref_path = COMP_DIR / "motif_mutation_crossref.csv"
    
    if not stats_path.exists() or not crossref_path.exists():
        print("[erreur] Fichiers requis introuvables. Vérifie les étapes précédentes.")
        return
        
    df_stats = pd.read_csv(stats_path)
    df_crossref = pd.read_csv(crossref_path)
    
    # 2. Fusion sur le nom du gène
    # On s'assure de nettoyer les noms de gènes si nécessaire
    merged_df = pd.merge(df_stats, df_crossref[['gene', 'warning_class']].drop_duplicates(), on='gene', how='inner')
    
    # 3. Agrégation par classe de toxicité
    toxicity_summary = merged_df.groupby('warning_class').agg(
        gene_count=('gene', 'count'),
        mean_alteration_pct=('alteration_frequency_pct', 'mean'),
        mean_abs_gistic=('mean_abs_value', 'mean')
    ).reset_index().sort_values(by='mean_alteration_pct', ascending=False)
    
    # 4. Sauvegarde
    out_path = CLINICAL_DIR / "toxicity_class_alteration_summary.csv"
    merged_df.to_csv(CLINICAL_DIR / "genes_with_toxicity_classes.csv", index=False)
    toxicity_summary.to_csv(out_path, index=False)
    
    print(f"\n[saved] Rapport par classe de toxicité enregistré : {out_path}")
    print("\n--- Synthèse par classe de toxicité ---")
    print(toxicity_summary.to_string(index=False))

if __name__ == "__main__":
    main()
