from pathlib import Path
import pandas as pd
import requests

PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"
RWR_DIR = PROC_DIR / "rwr_results"
CLINICAL_DIR = PROC_DIR / "cbioportal"
CLINICAL_DIR.mkdir(exist_ok=True, parents=True)
BASE_URL = "https://www.cbioportal.org/api"

def main():
    print("--- 1. Chargement des gènes cibles issus du RWR ---")
    rwr_path = RWR_DIR / "rwr_top_propagated_genes.csv"
    if not rwr_path.exists():
        print("[erreur] Le fichier RWR est introuvable. Lance d'abord 'python scripts/15_random_walk.py'.")
        return
        
    df_rwr = pd.read_csv(rwr_path)
    # Extraire les gènes cibles
    top_genes = df_rwr.groupby("warning_class").head(10)["gene"].unique().tolist()
    
    print("--- 2. Mapping des gènes ---")
    response = requests.get(f"{BASE_URL}/genes", params={"projection": "SUMMARY"})
    all_api_genes = response.json()
    gene_map = {g["hugoGeneSymbol"]: g.get("entrezGeneId") for g in all_api_genes if g.get("hugoGeneSymbol") in top_genes}
    
    print("--- 3. Récupération des cohortes PanCancer ---")
    studies_res = requests.get(f"{BASE_URL}/studies")
    all_studies = studies_res.json()
    pancan_studies = [s["studyId"] for s in all_studies if "pan_can_atlas_2018" in s["studyId"]]
    
    print(f"Extraction des données CNA pour {len(pancan_studies)} cohortes...")
    all_data = []

    for study_id in pancan_studies:
        cna_profile_id = f"{study_id}_gistic"
        sample_list_id = f"{study_id}_all"
        url = f"{BASE_URL}/molecular-profiles/{cna_profile_id}/molecular-data"
        
        for gene, entrez_id in gene_map.items():
            try:
                # Requête unitaire par gène pour éviter les erreurs 405/413
                res = requests.get(url, params={"entrezGeneId": entrez_id, "sampleListId": sample_list_id})
                if res.status_code == 200 and res.json():
                    for item in res.json():
                        all_data.append({
                            "gene": gene,
                            "study_id": study_id,
                            "sample_id": item.get("sampleId"),
                            "value": item.get("value")
                        })
            except Exception:
                continue

    if all_data:
        df_detailed = pd.DataFrame(all_data)
        out_path = CLINICAL_DIR / "pancancer_detailed_cna.csv"
        df_detailed.to_csv(out_path, index=False)
        
        # Résumé pour console
        print(f"\n[saved] Données détaillées exportées : {out_path}")
        print(f"Total des événements CNA enregistrés : {len(df_detailed)}")
        
        # Résumé rapide
        summary = df_detailed[df_detailed['value'].abs() >= 1].groupby('gene').size().reset_index(name='alt_count')
        print("\nTop 5 des gènes les plus altérés :")
        print(summary.sort_values(by='alt_count', ascending=False).head(5).to_string(index=False))
    else:
        print("[info] Aucune donnée CNA n'a pu être récupérée.")

if __name__ == "__main__":
    main()
