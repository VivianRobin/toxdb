from pathlib import Path
import pandas as pd

PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"
CLINICAL_DIR = PROC_DIR / "cbioportal"

def main():
    print("--- Génération du profil fonctionnel des toxicités ---")
    path = CLINICAL_DIR / "genes_with_toxicity_classes.csv"
    if not path.exists():
        print("[erreur] Fichier de synthèse introuvable.")
        return
        
    df = pd.read_csv(path)
    
    # Résumé par classe de toxicité et métriques moyennes
    summary = df.groupby('warning_class').agg(
        n_genes=('gene', 'count'),
        avg_freq_pct=('alteration_frequency_pct', 'mean'),
        top_gene=('gene', lambda x: x.iloc[0])
    ).reset_index()
    
    out_path = CLINICAL_DIR / "functional_toxicity_profile.csv"
    summary.to_csv(out_path, index=False)
    print(f"[saved] Profil fonctionnel enregistré : {out_path}")
    print(summary.to_string(index=False))

if __name__ == "__main__":
    main()
