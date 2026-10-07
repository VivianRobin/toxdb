"""
Ajoute un ou plusieurs médicaments PAR NOM (indépendamment du flag
warning_type='Withdrawn' de ChEMBL, qui n'est pas toujours renseigné
uniformément) et récupère leurs cibles :
  - mécanisme officiel (endpoint `mechanism`)
  - bioactivité off-target (endpoint `activity`, pchembl >= seuil)

Utile pour construire des cas de comparaison ciblés (ex: plusieurs
cardiotoxiques hERG connus, plusieurs hépatotoxiques...) sans dépendre
entièrement du filtre automatique de 02_download_chembl_toxic_drugs.py.

Usage :
  python scripts/02c_add_drug_by_name.py "Terfenadine" "Cisapride" "Sertindole"

Sortie : fusionne dans data/processed/toxic_drugs.csv et
         data/processed/drug_targets.csv / drug_targets_offtarget.csv
         (dédoublonnage automatique par chembl_id)
"""
import sys
from pathlib import Path

import pandas as pd
from chembl_webresource_client.new_client import new_client
from tqdm import tqdm

PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"
PCHEMBL_THRESHOLD = 6


def resolve_chembl_id(name: str) -> dict | None:
    molecule = new_client.molecule
    res = molecule.search(name)
    res = list(res)
    if not res:
        print(f"  [!] '{name}' introuvable dans ChEMBL")
        return None
    # on prend le meilleur match (le premier résultat de search, trié par pertinence)
    m = res[0]
    return {"chembl_id": m["molecule_chembl_id"], "pref_name": m.get("pref_name") or name.upper()}


def get_mechanism_targets(chembl_id: str) -> list[dict]:
    mechanism = new_client.mechanism
    target = new_client.target
    rows = []
    for m in mechanism.filter(molecule_chembl_id=chembl_id):
        tid = m.get("target_chembl_id")
        if not tid:
            continue
        t = list(target.filter(target_chembl_id=tid).only("target_chembl_id", "target_components"))
        gene_symbol = None
        if t and t[0].get("target_components"):
            comps = t[0]["target_components"]
            if comps:
                syns = comps[0].get("target_component_synonyms", [])
                gene_syns = [s["component_synonym"] for s in syns if s.get("syn_type") == "GENE_SYMBOL"]
                gene_symbol = gene_syns[0] if gene_syns else comps[0].get("accession")
        if gene_symbol:
            rows.append({"chembl_id": chembl_id, "target_gene_symbol": gene_symbol, "action_type": m.get("action_type")})
    return rows


def get_offtarget_activities(chembl_id: str) -> list[dict]:
    activity = new_client.activity
    target = new_client.target
    rows = []
    acts = activity.filter(
        molecule_chembl_id=chembl_id,
        pchembl_value__gte=PCHEMBL_THRESHOLD,
        assay_type="B",
        target_organism="Homo sapiens",
    ).only("target_chembl_id", "pchembl_value", "standard_type")

    for a in acts:
        tid = a.get("target_chembl_id")
        if not tid:
            continue
        t = list(target.filter(target_chembl_id=tid, target_type="SINGLE PROTEIN").only("target_components"))
        if not t:
            continue
        comps = t[0].get("target_components") or []
        gene_symbol = None
        if comps:
            syns = comps[0].get("target_component_synonyms", [])
            gene_syns = [s["component_synonym"] for s in syns if s.get("syn_type") == "GENE_SYMBOL"]
            gene_symbol = gene_syns[0] if gene_syns else comps[0].get("accession")
        if gene_symbol:
            rows.append(
                {
                    "chembl_id": chembl_id,
                    "target_gene_symbol": gene_symbol,
                    "pchembl_value": float(a["pchembl_value"]),
                    "standard_type": a.get("standard_type"),
                }
            )
    # garde la meilleure affinité par cible
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("pchembl_value", ascending=False).drop_duplicates(subset=["target_gene_symbol"])
    return df.to_dict("records")


def append_dedup(path: Path, new_rows: list[dict], subset: list[str]):
    new_df = pd.DataFrame(new_rows)
    if path.exists():
        old_df = pd.read_csv(path)
        combined = pd.concat([old_df, new_df], ignore_index=True).drop_duplicates(subset=subset)
    else:
        combined = new_df
    combined.to_csv(path, index=False)
    return combined


def main():
    if len(sys.argv) < 2:
        print('Usage: python scripts/02c_add_drug_by_name.py "Drug1" "Drug2" ...')
        sys.exit(1)

    names = sys.argv[1:]
    new_drugs, new_targets, new_offtargets = [], [], []

    for name in names:
        print(f"[{name}]")
        info = resolve_chembl_id(name)
        if not info:
            continue

        new_drugs.append(
            {
                "chembl_id": info["chembl_id"],
                "pref_name": info["pref_name"],
                "toxicity_reason": "manually_added",
                "warning_class": None,
                "warning_country": None,
                "warning_year": None,
            }
        )

        mech_targets = get_mechanism_targets(info["chembl_id"])
        print(f"  -> {len(mech_targets)} cible(s) mécanisme officiel")
        new_targets.extend(mech_targets)

        offtargets = get_offtarget_activities(info["chembl_id"])
        print(f"  -> {len(offtargets)} cible(s) off-target (bioactivité)")
        new_offtargets.extend(offtargets)

    append_dedup(PROC_DIR / "toxic_drugs.csv", new_drugs, subset=["chembl_id"])
    if new_targets:
        append_dedup(PROC_DIR / "drug_targets.csv", new_targets, subset=["chembl_id", "target_gene_symbol"])
    if new_offtargets:
        append_dedup(
            PROC_DIR / "drug_targets_offtarget.csv", new_offtargets, subset=["chembl_id", "target_gene_symbol"]
        )

    print("\n[done] CSVs mis à jour. Relance scripts/03_load_neo4j.py pour charger dans Neo4j.")


if __name__ == "__main__":
    main()
