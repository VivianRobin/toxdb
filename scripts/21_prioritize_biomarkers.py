from pathlib import Path
import pandas as pd

PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"
CLINICAL_DIR = PROC_DIR / "cbioportal"

def main():
    print("--- Calcul du score de priorisation des biomarqueurs ---")
    path = CLINICAL_DIR / "genes_with_toxicity_classes.csv"
    if not path.exists():
        return
        
    df = pd.read_csv(path)
    
    # Score composite : Normalisation min-max simple ou produit pondéré
    df['priority_score'] = (df['alteration_frequency_pct'] / 100) * df['mean_abs_value']
    
    # Sélection du Top 10 des biomarqueurs universels
    top_biomarkers = df.sort_values(by='priority_score', ascending=False).head(10)
    
    out_path = CLINICAL_DIR / "prioritized_biomarkers_panel.csv"
    top_biomarkers.to_csv(out_path, index=False)
    
    print(f"[saved] Panel de biomarqueurs priorisé : {out_path}")
    print(top_biomarkers[['gene', 'warning_class', 'alteration_frequency_pct', 'priority_score']].to_string(index=False))

if __name__ == "__main__":
    main()
