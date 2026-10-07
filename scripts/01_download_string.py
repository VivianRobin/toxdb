"""
Télécharge le réseau d'interactions protéine-protéine humain depuis STRING,
filtre par score de confiance, et sauvegarde une version "processed" prête
pour l'import Neo4j.

Source : https://string-db.org/cgi/download
Pas de compte / clé API nécessaire.

Sortie : data/processed/ppi_edges.csv
Colonnes : protein1, protein2, combined_score
(les IDs sont au format STRING, ex: 9606.ENSP00000000233 -- on les nettoiera
 à l'étape 05 lors du mapping vers des noms de gènes/UniProt)
"""
import gzip
import os
import shutil
from pathlib import Path

import pandas as pd
import requests
from dotenv import load_dotenv
from tqdm import tqdm

load_dotenv()

SPECIES = os.getenv("STRING_SPECIES", "9606")
VERSION = os.getenv("STRING_VERSION", "12.0")
SCORE_THRESHOLD = int(os.getenv("STRING_SCORE_THRESHOLD", "700"))

RAW_DIR = Path(__file__).resolve().parents[1] / "data" / "raw"
PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"
RAW_DIR.mkdir(parents=True, exist_ok=True)
PROC_DIR.mkdir(parents=True, exist_ok=True)

# Fichier "links" (scores agrégés) -- suffisant pour un réseau consensus.
# Le fichier "links.full" contient le détail par canal d'évidence (utile plus tard).
FILENAME = f"{SPECIES}.protein.links.v{VERSION}.txt.gz"
URL = f"https://stringdb-static.org/download/protein.links.v{VERSION}/{FILENAME}"

# Fichier d'info protéines : permet de récupérer le nom du gène (preferred_name)
INFO_FILENAME = f"{SPECIES}.protein.info.v{VERSION}.txt.gz"
INFO_URL = f"https://stringdb-static.org/download/protein.info.v{VERSION}/{INFO_FILENAME}"


def download(url: str, dest: Path):
    if dest.exists():
        print(f"[skip] {dest.name} déjà présent")
        return
    print(f"[download] {url}")
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        with open(dest, "wb") as f, tqdm(total=total, unit="B", unit_scale=True) as pbar:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
                pbar.update(len(chunk))


def main():
    links_gz = RAW_DIR / FILENAME
    info_gz = RAW_DIR / INFO_FILENAME

    download(URL, links_gz)
    download(INFO_URL, info_gz)

    print("[load] protein.info (mapping ID -> preferred_name)")
    with gzip.open(info_gz, "rt") as f:
        info = pd.read_csv(f, sep="\t")
    # colonnes attendues: #string_protein_id, preferred_name, protein_size, annotation
    info = info.rename(columns={"#string_protein_id": "string_protein_id"})
    id_to_name = dict(zip(info["string_protein_id"], info["preferred_name"]))

    print("[load] protein.links (edges)")
    with gzip.open(links_gz, "rt") as f:
        edges = pd.read_csv(f, sep=" ")

    n_before = len(edges)
    edges = edges[edges["combined_score"] >= SCORE_THRESHOLD].copy()
    print(f"[filter] score >= {SCORE_THRESHOLD}: {n_before} -> {len(edges)} interactions")

    edges["gene1"] = edges["protein1"].map(id_to_name)
    edges["gene2"] = edges["protein2"].map(id_to_name)
    edges = edges.dropna(subset=["gene1", "gene2"])

    out = edges[["gene1", "gene2", "combined_score"]].rename(
        columns={"gene1": "source", "gene2": "target", "combined_score": "score"}
    )
    out_path = PROC_DIR / "ppi_edges.csv"
    out.to_csv(out_path, index=False)
    print(f"[done] {len(out)} arêtes -> {out_path}")


if __name__ == "__main__":
    main()
