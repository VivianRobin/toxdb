from pathlib import Path
import pandas as pd

PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"
CLINICAL_DIR = PROC_DIR / "cbioportal"

def main():
    print("--- Analyse statistique des altérations CNA ---")
    data_path = CLINICAL_DIR / "pancancer_detailed_cna.csv"
    if not data_path.exists():
        print("[erreur] Le fichier pancancer_detailed_cna.csv est introuvable.")
        return
        
    df = pd.read_csv(data_path)
    
    # Définition d'une altération significative (Gain >= 1 ou Perte <= -1)
    df['is_altered'] = df['value'].abs() >= 1
    
    # Calcul des métriques par gène
    stats = df.groupby('gene').agg(
        total_samples=('sample_id', 'count'),
        altered_samples=('is_altered', 'sum'),
        mean_abs_value=('value', lambda x: x.abs().mean()),
        mean_value=('value', 'mean')
    ).reset_index()
    
    # Calcul du pourcentage d'altération
    stats['alteration_frequency_pct'] = (stats['altered_samples'] / stats['total_samples']) * 100
    
    # Tri par fréquence d'altération décroissante
    stats_sorted = stats.sort_values(by='alteration_frequency_pct', ascending=False)
    
    # Sauvegarde des résultats
    out_csv = CLINICAL_DIR / "top_altered_genes_summary.csv"
    stats_sorted.to_csv(out_csv, index=False)
    
    print(f"\n[saved] Rapport statistique enregistré : {out_csv}")
    print("\n--- Top 10 des gènes les plus altérés (Fréquence PanCancer) ---")
    print(stats_sorted[['gene', 'alteration_frequency_pct', 'mean_abs_value']].head(10).to_string(index=False))

if __name__ == "__main__":
    main()
