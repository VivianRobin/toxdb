"""
Récupère depuis ChEMBL :
  1. La liste des médicaments approuvés puis RETIRÉS du marché pour raisons
     de toxicité (endpoint drug_warning, champ warning_class).
  2. Pour chacun, les cibles protéiques connues (via mechanism -> target).

Aucune clé API / compte nécessaire (web service public EBI).

Sorties :
  data/processed/toxic_drugs.csv        (chembl_id, pref_name, warning_class, warning_country, warning_year)
  data/processed/drug_targets.csv       (chembl_id, drug_name, target_chembl_id, target_gene_symbol, action_type)
"""
from pathlib import Path

import pandas as pd
from chembl_webresource_client.new_client import new_client
from tqdm import tqdm

PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"
PROC_DIR.mkdir(parents=True, exist_ok=True)


def get_withdrawn_drugs() -> pd.DataFrame:
    drug_warning = new_client.drug_warning
    molecule = new_client.molecule

    # warning_type == 'Withdrawn' isole les retraits (par opposition aux simples
    # boxed warnings qui laissent le médicament sur le marché)
    results = drug_warning.filter(warning_type="Withdrawn")
    rows = []
    for r in tqdm(results, desc="drug_warning (withdrawn)"):
        rows.append(
            {
                "chembl_id": r.get("parent_molecule_chembl_id"),
                # efo_term = motif/type de toxicité (ex: "hepatotoxicity"), PAS
                # le nom du médicament -- on le garde à part, distinct de pref_name.
                "toxicity_reason": r.get("efo_term"),
                "warning_class": r.get("warning_class"),
                "warning_country": r.get("warning_country"),
                "warning_year": r.get("warning_year"),
            }
        )
    df = pd.DataFrame(rows).drop_duplicates(subset=["chembl_id"])

    # Le vrai nom du médicament se récupère via l'endpoint molecule, pas
    # drug_warning. On le résout ici, en cache par chembl_id.
    print("Résolution des noms de médicaments via l'endpoint molecule...")
    name_cache = {}
    for cid in tqdm(df["chembl_id"].dropna().unique(), desc="molecule -> pref_name"):
        m = molecule.filter(molecule_chembl_id=cid).only("molecule_chembl_id", "pref_name")
        m = list(m)
        name_cache[cid] = m[0]["pref_name"] if m and m[0].get("pref_name") else cid

    df["pref_name"] = df["chembl_id"].map(name_cache)
    return df


def get_targets_for_drugs(chembl_ids: list[str]) -> pd.DataFrame:
    mechanism = new_client.mechanism
    target = new_client.target

    rows = []
    target_cache = {}

    for cid in tqdm(chembl_ids, desc="mechanism -> target"):
        if not cid:
            continue
        mechs = mechanism.filter(molecule_chembl_id=cid)
        for m in mechs:
            tid = m.get("target_chembl_id")
            if not tid:
                continue
            if tid not in target_cache:
                t = target.filter(target_chembl_id=tid).only(
                    "target_chembl_id", "pref_name", "target_components"
                )
                t = list(t)
                target_cache[tid] = t[0] if t else None
            t_info = target_cache[tid]
            gene_symbol = None
            if t_info and t_info.get("target_components"):
                comps = t_info["target_components"]
                if comps:
                    # target_component_synonyms peut contenir le gene symbol
                    syns = comps[0].get("target_component_synonyms", [])
                    gene_syns = [s["component_synonym"] for s in syns if s.get("syn_type") == "GENE_SYMBOL"]
                    gene_symbol = gene_syns[0] if gene_syns else comps[0].get("accession")

            rows.append(
                {
                    "chembl_id": cid,
                    "target_chembl_id": tid,
                    "target_name": t_info.get("pref_name") if t_info else None,
                    "target_gene_symbol": gene_symbol,
                    "action_type": m.get("action_type"),
                }
            )
    return pd.DataFrame(rows)


def main():
    print("[1/2] Récupération des médicaments retirés pour toxicité...")
    toxic_drugs = get_withdrawn_drugs()
    toxic_drugs.to_csv(PROC_DIR / "toxic_drugs.csv", index=False)
    print(f"  -> {len(toxic_drugs)} médicaments -> data/processed/toxic_drugs.csv")

    print("[2/2] Récupération des cibles protéiques associées...")
    targets = get_targets_for_drugs(toxic_drugs["chembl_id"].dropna().unique().tolist())
    targets = targets.dropna(subset=["target_gene_symbol"])
    targets.to_csv(PROC_DIR / "drug_targets.csv", index=False)
    print(f"  -> {len(targets)} liens médicament-cible -> data/processed/drug_targets.csv")


if __name__ == "__main__":
    main()
