"""
Récupère la fréquence de mutation pan-cancer (TCGA PanCancer Atlas, exome
complet) pour chaque protéine du réseau IPP, via l'API REST publique de
cBioPortal (aucune clé/compte nécessaire).

Pourquoi TCGA PanCancer Atlas plutôt que MSK-IMPACT : MSK-IMPACT est un
PANEL CIBLÉ (~400 gènes) -- un gène absent du panel apparaîtrait à tort
comme "jamais muté" alors qu'il n'a simplement jamais été testé. Les 32
études PanCancer Atlas sont en exome complet, donc chaque gène du génome
est uniformément couvert -- essentiel pour ne pas biaiser la comparaison.

Les IDs des 32 études sont découverts DYNAMIQUEMENT via /api/studies
(pattern 'xxx_tcga_pan_can_atlas_2018'), pas codés en dur -- reproductible
même si la liste change côté cBioPortal.

Étapes :
  1. Découverte des études PanCancer Atlas + taille de cohorte (échantillons)
  2. Résolution symbole de gène -> Entrez Gene ID (endpoint /genes/fetch, en lot)
  3. Pour chaque étude, récupération des mutations pour nos gènes (paginé)
  4. Agrégation : fréquence = (# échantillons mutés, toutes études) / (# total échantillons)

Usage :
  # Sur toutes les protéines déjà chargées dans Neo4j (recommandé)
  python scripts/08_download_cbioportal_mutations.py --from-neo4j

  # Sur une liste explicite (test rapide)
  python scripts/08_download_cbioportal_mutations.py --genes KCNH2,KCNQ1,SCN5A,TP53

Sortie : data/processed/gene_mutation_frequency.csv
Colonnes : gene, entrez_id, n_mutated_samples, n_total_samples, mutation_frequency
"""
import argparse
import os
import time
from pathlib import Path

import pandas as pd
import requests
from dotenv import load_dotenv
from tqdm import tqdm

load_dotenv()

BASE_URL = "https://www.cbioportal.org/api"
PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"

GENE_CHUNK_SIZE = 400          # taille des lots pour la résolution de gènes / requêtes mutations
PAGE_SIZE = 10000              # pagination des résultats de mutations
MAX_PAGES_PER_CHUNK = 200      # garde-fou anti-boucle infinie
REQUEST_TIMEOUT = 60
RETRY_ATTEMPTS = 3


def _post_with_retry(url: str, **kwargs) -> requests.Response:
    for attempt in range(RETRY_ATTEMPTS):
        try:
            r = requests.post(url, timeout=REQUEST_TIMEOUT, **kwargs)
            r.raise_for_status()
            return r
        except requests.RequestException as e:
            if attempt == RETRY_ATTEMPTS - 1:
                body = getattr(e.response, "text", "")[:500] if getattr(e, "response", None) is not None else ""
                print(f"  [error body] {body}")
                raise
            print(f"  [retry {attempt + 1}/{RETRY_ATTEMPTS}] {e}")
            time.sleep(2 * (attempt + 1))


def discover_pancan_atlas_studies() -> list[dict]:
    """Découvre les 32 études PanCancer Atlas dynamiquement (pas de liste codée en dur)."""
    print("[1/4] Découverte des études TCGA PanCancer Atlas...")
    r = requests.get(f"{BASE_URL}/studies", params={"projection": "SUMMARY"}, timeout=REQUEST_TIMEOUT)
    r.raise_for_status()
    all_studies = r.json()

    pancan_studies = [s for s in all_studies if s["studyId"].endswith("_tcga_pan_can_atlas_2018")]
    print(f"  -> {len(pancan_studies)} études trouvées (attendu: ~32)")
    if len(pancan_studies) == 0:
        raise RuntimeError(
            "Aucune étude PanCancer Atlas trouvée -- le pattern d'ID a peut-être changé côté cBioPortal. "
            "Vérifie manuellement sur https://www.cbioportal.org/datasets"
        )
    return pancan_studies


def get_sample_list_size(study_id: str) -> int:
    sample_list_id = f"{study_id}_all"
    r = requests.get(f"{BASE_URL}/sample-lists/{sample_list_id}", timeout=REQUEST_TIMEOUT)
    if r.status_code != 200:
        return 0
    return len(r.json().get("sampleIds", []))


def resolve_entrez_ids(gene_symbols: list[str]) -> dict:
    """Résout HUGO symbol -> entrezGeneId, en lot."""
    print(f"[2/4] Résolution de {len(gene_symbols)} symboles de gènes -> Entrez ID...")
    mapping = {}
    unique_symbols = list(dict.fromkeys(gene_symbols))  # dédoublonne en gardant l'ordre

    for i in tqdm(range(0, len(unique_symbols), GENE_CHUNK_SIZE), desc="genes/fetch"):
        chunk = unique_symbols[i : i + GENE_CHUNK_SIZE]
        r = _post_with_retry(
            f"{BASE_URL}/genes/fetch",
            params={"geneIdType": "HUGO_GENE_SYMBOL"},
            json=chunk,
            headers={"Content-Type": "application/json"},
        )
        for gene in r.json():
            mapping[gene["hugoGeneSymbol"]] = gene["entrezGeneId"]

    print(f"  -> {len(mapping)}/{len(unique_symbols)} symboles résolus")
    return mapping


def fetch_mutated_samples_for_study(study_id: str, entrez_ids: list[int]) -> dict:
    """Retourne {entrezGeneId: set(sampleId muté)} pour une étude donnée."""
    profile_id = f"{study_id}_mutations"
    sample_list_id = f"{study_id}_all"
    result = {eid: set() for eid in entrez_ids}

    for i in range(0, len(entrez_ids), GENE_CHUNK_SIZE):
        gene_chunk = entrez_ids[i : i + GENE_CHUNK_SIZE]
        page = 0
        while page < MAX_PAGES_PER_CHUNK:
            r = _post_with_retry(
                f"{BASE_URL}/molecular-profiles/{profile_id}/mutations/fetch",
                params={
                    "projection": "ID",  # champs minimaux, réponse plus légère/rapide
                    "pageSize": PAGE_SIZE,
                    "pageNumber": page,
                },
                json={"entrezGeneIds": gene_chunk, "sampleListId": sample_list_id},
                headers={"Content-Type": "application/json"},
            )
            data = r.json()
            if not data:
                break
            for mut in data:
                eid = mut.get("entrezGeneId")
                sid = mut.get("sampleId")
                if eid in result and sid:
                    result[eid].add(sid)
            if len(data) < PAGE_SIZE:
                break
            page += 1

    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--from-neo4j", action="store_true", help="Utilise toutes les protéines déjà en base")
    src.add_argument("--genes", type=str, help="Liste de gènes séparés par virgules (test rapide)")
    args = parser.parse_args()

    if args.genes:
        gene_symbols = [g.strip() for g in args.genes.split(",") if g.strip()]
    else:
        from neo4j import GraphDatabase

        uri = os.getenv("NEO4J_URI", "bolt://localhost:7687")
        user = os.getenv("NEO4J_USER", "neo4j")
        password = os.getenv("NEO4J_PASSWORD", "toxkg-dev-pw")
        driver = GraphDatabase.driver(uri, auth=(user, password))
        with driver.session() as session:
            gene_symbols = [r["name"] for r in session.run("MATCH (p:Protein) RETURN p.name AS name")]
        driver.close()
        print(f"[from-neo4j] {len(gene_symbols)} protéines récupérées depuis Neo4j")

    studies = discover_pancan_atlas_studies()

    print("[setup] taille des cohortes par étude...")
    total_samples = 0
    for s in tqdm(studies, desc="sample-lists"):
        n = get_sample_list_size(s["studyId"])
        s["n_samples"] = n
        total_samples += n
    print(f"  -> {total_samples} échantillons au total (pan-cancer)")

    symbol_to_entrez = resolve_entrez_ids(gene_symbols)
    entrez_to_symbol = {v: k for k, v in symbol_to_entrez.items()}
    entrez_ids = list(entrez_to_symbol.keys())

    print(f"[3/4] Récupération des mutations sur {len(studies)} études...")
    mutated_samples = {eid: set() for eid in entrez_ids}  # union sur toutes études (studyId::sampleId)

    for s in tqdm(studies, desc="studies"):
        study_result = fetch_mutated_samples_for_study(s["studyId"], entrez_ids)
        for eid, sample_ids in study_result.items():
            # préfixe par studyId pour éviter les collisions d'ID d'échantillon entre études
            mutated_samples[eid].update(f"{s['studyId']}::{sid}" for sid in sample_ids)

    print("[4/4] Calcul des fréquences...")
    rows = []
    for eid in entrez_ids:
        n_mut = len(mutated_samples[eid])
        rows.append(
            {
                "gene": entrez_to_symbol[eid],
                "entrez_id": eid,
                "n_mutated_samples": n_mut,
                "n_total_samples": total_samples,
                "mutation_frequency": n_mut / total_samples if total_samples else 0.0,
            }
        )

    df = pd.DataFrame(rows).sort_values("mutation_frequency", ascending=False)
    PROC_DIR.mkdir(parents=True, exist_ok=True)
    out_path = PROC_DIR / "gene_mutation_frequency.csv"
    df.to_csv(out_path, index=False)
    print(f"[done] {len(df)} gènes -> {out_path}")
    print(df.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
