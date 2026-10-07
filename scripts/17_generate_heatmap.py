from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Force le mode sans affichage graphique (anti-erreur TclError)
import matplotlib.pyplot as plt
import seaborn as sns

PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"
CLINICAL_DIR = PROC_DIR / "cbioportal"

def main():
    print("--- Chargement des données détaillées CNA ---")
    data_path = CLINICAL_DIR / "pancancer_detailed_cna.csv"
    if not data_path.exists():
        print("[erreur] Le fichier pancancer_detailed_cna.csv est introuvable. Relance le script 16.")
        return
        
    df = pd.read_csv(data_path)
    
    # Création de la matrice Gène x Étude (moyenne des valeurs GISTIC)
    mat = df.pivot_table(index="gene", columns="study_id", values="value", aggfunc="mean").fillna(0)
    
    print(f"Matrice prête : {mat.shape[0]} gènes x {mat.shape[1]} cohortes.")
    
    # Génération de la Heatmap
    plt.figure(figsize=(14, 10))
    sns.heatmap(mat, cmap="vlag", center=0, cbar_kws={'label': 'Mean GISTIC CNA Value'},
                linewidths=0.5, linecolor='gray', 
                xticklabels=True, yticklabels=True)
    
    plt.title("Altérations CNA (GISTIC) - TCGA PanCancer Atlas", fontsize=14, fontweight='bold')
    plt.xlabel("Cohortes TCGA", fontsize=12)
    plt.ylabel("Gènes de Toxicité", fontsize=12)
    plt.xticks(rotation=90, fontsize=8)
    plt.yticks(fontsize=9)
    
    plt.tight_layout()
    out_pdf = CLINICAL_DIR / "toxicity_cna_heatmap.pdf"
    plt.savefig(out_pdf, dpi=300)
    print(f"[saved] Heatmap enregistrée avec succès dans : {out_pdf}")

if __name__ == "__main__":
    main()
