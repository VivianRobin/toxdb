from pathlib import Path
import pandas as pd
import gseapy as gp

PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"
COMPARISON_DIR = PROC_DIR / "comparisons"
ENRICH_DIR = PROC_DIR / "enrichments"
ENRICH_DIR.mkdir(exist_ok=True, parents=True)

def main():
    # 1. Charger tes gènes priorisés
    crossref = pd.read_csv(COMPARISON_DIR / "motif_mutation_crossref.csv")
    
    classes = crossref["warning_class"].unique()
    
    for cls in classes:
        print(f"\n--- Analyse d'enrichissement pour : {cls} ---")
        
        # Isoler le Top 20 des gènes de la classe
        df_cls = crossref[crossref["warning_class"] == cls].head(20)
        gene_list = df_cls["gene"].tolist()
        
        if not gene_list:
            continue

        # 2. Interroger les bases de signatures (ex: MSigDB Hallmark & KEGG)
        # gseapy permet d'interroger Enrichr directement en ligne
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=['MSigDB_Hallmark_2020', 'KEGG_2021_Human'],
            organism='human',
            outdir=None # Mode programmatique
        )
        
        res_df = enr.res2d
        
        # Filtrer les résultats significatifs (p-value ajustée < 0.05)
        sig_res = res_df[res_df['Adjusted P-value'] < 0.05]
        
        if not sig_res.empty:
            out_path = ENRICH_DIR / f"enrichment_{cls.replace(' ', '_')}.csv"
            sig_res.to_csv(out_path, index=False)
            print(f"[saved] {out_path} ({len(sig_res)} signatures significatives trouvées)")
            
            # Afficher le top 3 des signatures les plus touchées
            print(sig_res[['Gene_set', 'Term', 'Adjusted P-value']].head(3).to_string(index=False))
        else:
            print("[info] Aucune signature hautement significative trouvée pour ce seuil.")

if __name__ == "__main__":
    main()
