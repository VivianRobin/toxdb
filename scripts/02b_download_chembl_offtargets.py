"""
Complète drug_targets.csv avec les cibles OFF-TARGET des médicaments toxiques,
via l'endpoint `activity` de ChEMBL (bioactivité mesurée), en plus des cibles
"officielles" déjà récupérées par 02_download_chembl_toxic_drugs.py (endpoint
`mechanism`, qui ne donne que la cible pharmacologique principale).

Pourquoi : la cause de toxicité est souvent une liaison SECONDAIRE non
recherchée (ex: blocage du canal hERG/KCNH2 par l'astémizole, un
antihistaminique dont la cible primaire est HRH1). Le mécanisme "officiel"
ne capture pas ça -- il faut regarder toute la bioactivité mesurée.

Filtre appliqué :
  - pchembl_value >= PCHEMBL_THRESHOLD (6 = activité submicromolaire, IC50/Ki/Kd < 1µM)
  - assay_type == 'B' (binding assay, le plus direct pour une cible protéique)
  - target_type == 'SINGLE PROTEIN' (on exclut cell lines, organismes entiers...)

Aucune clé API nécessaire.

Sortie : data/processed/drug_targets_offtarget.csv
Colonnes : chembl_id, target_chembl_id, target_gene_symbol, pchembl_value, standard_type
"""
from pathlib import Path

import pandas as pd
from chembl_webresource_client.new_client import new_client
from tqdm import tqdm

PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"

PCHEMBL_THRESHOLD = 6  # correspond à IC50/Ki/Kd < 1 µM, seuil classique "activité significative"


def get_offtarget_activities(chembl_ids: list[str]) -> pd.DataFrame:
    activity = new_client.activity
    target = new_client.target
    target_cache = {}

    rows = []
    for cid in tqdm(chembl_ids, desc="activity (off-target bioactivity)"):
        if not cid:
            continue
        acts = activity.filter(
            molecule_chembl_id=cid,
            pchembl_value__gte=PCHEMBL_THRESHOLD,
            assay_type="B",
            target_organism="Homo sapiens",
        ).only(
            "target_chembl_id", "pchembl_value", "standard_type", "target_organism"
        )

        for a in acts:
            tid = a.get("target_chembl_id")
            if not tid:
                continue

            if tid not in target_cache:
                t = target.filter(target_chembl_id=tid, target_type="SINGLE PROTEIN").only(
                    "target_chembl_id", "pref_name", "target_components"
                )
                t = list(t)
                target_cache[tid] = t[0] if t else None

            t_info = target_cache[tid]
            if t_info is None:
                continue  # pas une SINGLE PROTEIN -> on ignore (cell line, organisme entier, etc.)

            gene_symbol = None
            comps = t_info.get("target_components") or []
            if comps:
                syns = comps[0].get("target_component_synonyms", [])
                gene_syns = [s["component_synonym"] for s in syns if s.get("syn_type") == "GENE_SYMBOL"]
                gene_symbol = gene_syns[0] if gene_syns else comps[0].get("accession")

            if not gene_symbol:
                continue

            rows.append(
                {
                    "chembl_id": cid,
                    "target_chembl_id": tid,
                    "target_gene_symbol": gene_symbol,
                    "pchembl_value": a.get("pchembl_value"),
                    "standard_type": a.get("standard_type"),
                }
            )

    df = pd.DataFrame(rows)
    if not df.empty:
        # une même paire (drug, target) peut avoir plusieurs mesures -> on garde la plus forte affinité
        df["pchembl_value"] = df["pchembl_value"].astype(float)
        df = df.sort_values("pchembl_value", ascending=False).drop_duplicates(
            subset=["chembl_id", "target_gene_symbol"]
        )
    return df


def main():
    toxic_drugs = pd.read_csv(PROC_DIR / "toxic_drugs.csv")
    chembl_ids = toxic_drugs["chembl_id"].dropna().unique().tolist()

    print(f"Recherche de bioactivités off-target pour {len(chembl_ids)} médicaments...")
    offtargets = get_offtarget_activities(chembl_ids)

    out_path = PROC_DIR / "drug_targets_offtarget.csv"
    offtargets.to_csv(out_path, index=False)
    print(f"[done] {len(offtargets)} liens médicament-cible (off-target) -> {out_path}")


if __name__ == "__main__":
    main()
