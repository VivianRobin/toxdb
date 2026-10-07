"""
Ré-analyse de la classe "hepatotoxicity" avec des listes de médicaments
CURÉES (externes à ChEMBL), pour tester si le signal non significatif
obtenu en 3.3/Table 2 (hepatotoxicity : n=32 drogues via le flag ChEMBL
warning_type='Withdrawn' + warning_class, aucun partage de motif au-delà
du hasard) est un vrai résultat biologique ou un artefact du label bruité.

Contexte (limitation citée dans le manuscrit, section Discussion) :
  "the automatic withdrawal flag is not uniformly populated in ChEMBL" /
  "the ChEMBL warning-class taxonomy is coarse and unevenly populated"
-> on remplace ce label par deux listes hépatotoxicité curées indépendamment
   (list_Drug_Hepatotoxicity.txt, 115 noms ; list_hepa.txt, 89 noms ;
   intersection = 72 noms "core" haute confiance, union = 132 "broad"
   haut rappel), et on refait exactement le même pipeline que
   10_motif_analysis.py (PageRank personnalisé -> sous-graphe induit ->
   cliques maximales -> motif nodes -> Jaccard -> permutation d'étiquettes),
   mais restreint à ce sous-ensemble curé plutôt qu'au warning_class ChEMBL.

Ce script NE MODIFIE PAS warning_class ni le flag d.toxic existants (pour
ne pas perturber 10/11/13/19 qui en dépendent) : il ajoute une propriété
séparée d.hepatox_curated ('core' | 'broad') sur les nœuds Drug concernés,
et fetch les cibles ChEMBL (mécanisme + off-target) pour les médicaments
curés qui ne seraient pas déjà dans le graphe (même logique que
02c_add_drug_by_name.py).

Étapes :
  1. Parsing des 2 listes -> noms normalisés, core = intersection, broad = union
  2. Résolution ChEMBL (molecule.search) de chaque nom -> chembl_id, avec
     cache CSV (reprise possible) et export d'une liste de noms NON résolus
     pour revue manuelle (certaines entrées des listes sources sont des
     codes de recherche ou coquilles, ex. "A77 1726", "NSC 352122",
     "cyclosporineCYCLOSPORINE")
  3. Pour chaque chembl_id résolu : si le Drug n'a pas encore de relation
     TARGETS dans Neo4j, fetch mécanisme + off-target (pChEMBL >= seuil) et
     charge-les ; dans tous les cas, tag d.hepatox_curated
  4. Motifs (PageRank + cliques), réutilisés depuis drug_motif_nodes.csv
     si déjà calculés par 10_motif_analysis.py, sinon calculés ici
  5. Test par permutation d'étiquettes (identique en logique à 10, mais
     un seul groupe testé à la fois : "core" et "broad" séparément) contre
     un fond composé de TOUS les médicaments toxiques ayant des cibles
     connues (curés ou non)
  6. Comparaison imprimée avec le résultat "hepatotoxicity" original
     (motif_enrichment_by_class.csv, label ChEMBL brut) si disponible

Usage :
  python scripts/22_hepatox_curated_reanalysis.py
  python scripts/22_hepatox_curated_reanalysis.py \
      --list1 data/raw/hepatox_curated/list_Drug_Hepatotoxicity.txt \
      --list2 data/raw/hepatox_curated/list_hepa.txt \
      --top-n 40 --min-clique-size 3 --n-permutations 1000

Sorties (data/processed/comparisons/) :
  - hepatox_curated_chembl_map.csv        nom -> chembl_id (+ statut résolution)
  - hepatox_curated_unresolved.csv        noms non résolus, à curer à la main
  - hepatox_curated_drug_motif_nodes.csv motif nodes des médicaments curés
                                           nouvellement calculés
  - hepatox_curated_motif_enrichment.csv résultat core vs broad vs (rappel)
                                           hepatotoxicity ChEMBL brut
"""
import argparse
import os
import time
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from dotenv import load_dotenv
from neo4j import GraphDatabase
from tqdm import tqdm

load_dotenv()

URI = os.getenv("NEO4J_URI", "bolt://localhost:7687")
USER = os.getenv("NEO4J_USER", "neo4j")
PASSWORD = os.getenv("NEO4J_PASSWORD", "toxkg-dev-pw")

RAW_DIR = Path(__file__).resolve().parents[1] / "data" / "raw" / "hepatox_curated"
PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"
COMPARISON_DIR = PROC_DIR / "comparisons"

GRAPH_NAME = "ppi"
PCHEMBL_THRESHOLD = 6
CHEMBL_MAX_RETRIES = 5
CHEMBL_RETRY_BACKOFF_S = 10  # doublé à chaque essai : 10, 20, 40, 80, 160s


def _get_chembl_client():
    """Import paresseux + retry avec backoff exponentiel : l'API ChEMBL/EBI
    renvoie régulièrement des 500 transitoires (incident serveur, pas un
    problème de requête). On retente au lieu de planter tout le script."""
    last_err = None
    for attempt in range(1, CHEMBL_MAX_RETRIES + 1):
        try:
            from chembl_webresource_client.new_client import new_client
            return new_client
        except Exception as e:  # noqa: BLE001
            last_err = e
            wait = CHEMBL_RETRY_BACKOFF_S * attempt
            print(f"  [!] ChEMBL indisponible (essai {attempt}/{CHEMBL_MAX_RETRIES}: {e.__class__.__name__}) -> nouvelle tentative dans {wait}s")
            time.sleep(wait)
    raise RuntimeError(
        f"ChEMBL API inaccessible après {CHEMBL_MAX_RETRIES} tentatives (dernier incident probablement côté EBI, "
        f"vérifie https://www.ebi.ac.uk/chembl/ dans un navigateur). Erreur : {last_err}"
    )


# --------------------------------------------------------------------------
# 1. Chargement + normalisation des listes curées
# --------------------------------------------------------------------------
def load_name_list(path: Path) -> set[str]:
    names = set()
    for line in path.read_text().splitlines():
        n = line.strip().lower()
        if not n:
            continue
        n = n.split(",")[0].strip()
        names.add(n)
    return names


def build_curated_sets(list1: Path, list2: Path) -> pd.DataFrame:
    a = load_name_list(list1)
    b = load_name_list(list2)
    core = a & b
    broad = a | b
    print(f"[curated] {list1.name}: {len(a)} noms | {list2.name}: {len(b)} noms")
    print(f"[curated] core (intersection) = {len(core)} | broad (union) = {len(broad)}")

    rows = [{"name": n, "hepatox_curated": "core"} for n in sorted(core)]
    rows += [{"name": n, "hepatox_curated": "broad"} for n in sorted(broad - core)]
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# 2. Résolution ChEMBL (avec cache, reprise possible)
# --------------------------------------------------------------------------
def resolve_chembl_id(name: str) -> dict | None:
    molecule = _get_chembl_client().molecule
    res = list(molecule.search(name))
    if not res:
        return None
    m = res[0]
    return {"chembl_id": m["molecule_chembl_id"], "pref_name": m.get("pref_name") or name.upper()}


def resolve_all(curated_df: pd.DataFrame, cache_path: Path) -> pd.DataFrame:
    cached = pd.read_csv(cache_path) if cache_path.exists() else pd.DataFrame(columns=["name", "chembl_id", "pref_name", "resolved"])
    already_resolved = set(cached.loc[cached["resolved"] == True, "name"])
    cached = cached[~cached["name"].isin(cached.loc[cached["resolved"] != True, "name"])]
    todo = curated_df[~curated_df["name"].isin(already_resolved)]

    new_rows = []
    for name in tqdm(todo["name"], desc="résolution ChEMBL"):
        try:
            info = resolve_chembl_id(name)
        except Exception as e:  # noqa: BLE001
            print(f"  [!] échec pour '{name}' ({e.__class__.__name__}: {e}) -> marqué non résolu, relance le script pour réessayer")
            info = None
        if info:
            new_rows.append({"name": name, "chembl_id": info["chembl_id"], "pref_name": info["pref_name"], "resolved": True})
        else:
            new_rows.append({"name": name, "chembl_id": None, "pref_name": None, "resolved": False})
        
        pd.concat([cached, pd.DataFrame(new_rows)], ignore_index=True).drop_duplicates(subset=["name"]).to_csv(cache_path, index=False)

    merged = pd.concat([cached, pd.DataFrame(new_rows)], ignore_index=True).drop_duplicates(subset=["name"])
    
    # Évite le conflit de noms de colonnes (_x, _y) lors du merge
    if "hepatox_curated" in merged.columns:
        merged = merged.drop(columns=["hepatox_curated"])

    merged = merged.merge(curated_df, on="name", how="right")
    merged.to_csv(cache_path, index=False)

    n_ok = merged["resolved"].sum()
    print(f"[resolve] {n_ok}/{len(merged)} noms résolus dans ChEMBL")
    return merged


# --------------------------------------------------------------------------
# 3. Cibles ChEMBL (mécanisme + off-target) pour les médicaments manquants
# --------------------------------------------------------------------------
def drug_has_targets(session, chembl_id: str) -> bool:
    r = session.run(
        "MATCH (:Drug {chembl_id: $id})-[:TARGETS]->() RETURN count(*) AS n", id=chembl_id
    ).single()
    return r["n"] > 0


def get_mechanism_targets(chembl_id: str) -> list[dict]:
    client = _get_chembl_client()
    mechanism = client.mechanism
    target = client.target
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


def get_offtarget_activities(chembl_id: str, pchembl_threshold: int) -> list[dict]:
    client = _get_chembl_client()
    activity = client.activity
    target = client.target
    rows = []
    acts = activity.filter(
        molecule_chembl_id=chembl_id,
        pchembl_value__gte=pchembl_threshold,
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
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("pchembl_value", ascending=False).drop_duplicates(subset=["target_gene_symbol"])
    return df.to_dict("records")


def fetch_and_load_missing_targets(session, chembl_ids: list[str], pchembl_threshold: int):
    missing = [cid for cid in chembl_ids if not drug_has_targets(session, cid)]
    print(f"[targets] {len(missing)}/{len(chembl_ids)} médicaments curés sans cibles connues -> fetch ChEMBL")

    for chembl_id in tqdm(missing, desc="fetch cibles manquantes"):
        mech = get_mechanism_targets(chembl_id)
        off = get_offtarget_activities(chembl_id, pchembl_threshold)
        for row in mech:
            session.run(
                """
                MERGE (d:Drug {chembl_id: $chembl_id})
                MERGE (p:Protein {name: $gene})
                MERGE (d)-[r:TARGETS]->(p)
                SET r.action_type = $action_type, r.source = 'mechanism'
                """,
                chembl_id=chembl_id, gene=row["target_gene_symbol"], action_type=row.get("action_type"),
            )
        for row in off:
            session.run(
                """
                MERGE (d:Drug {chembl_id: $chembl_id})
                MERGE (p:Protein {name: $gene})
                MERGE (d)-[r:TARGETS]->(p)
                SET r.source = coalesce(r.source, 'bioactivity'),
                    r.pchembl_value = $pchembl_value, r.standard_type = $standard_type
                """,
                chembl_id=chembl_id, gene=row["target_gene_symbol"],
                pchembl_value=row["pchembl_value"], standard_type=row.get("standard_type"),
            )


def tag_curated_drugs(session, resolved_df: pd.DataFrame):
    rows = resolved_df[resolved_df["resolved"] == True][["chembl_id", "pref_name", "hepatox_curated"]].to_dict("records")
    session.run(
        """
        UNWIND $rows AS row
        MERGE (d:Drug {chembl_id: row.chembl_id})
        ON CREATE SET d.name = row.pref_name, d.toxic = true, d.toxicity_reason = 'hepatox_curated_external'
        SET d.hepatox_curated = row.hepatox_curated
        """,
        rows=rows,
    )
    print(f"[tag] {len(rows)} médicaments tagués d.hepatox_curated dans Neo4j")


# --------------------------------------------------------------------------
# 4. Motifs (identique à 10_motif_analysis.py) + réutilisation du cache
# --------------------------------------------------------------------------
def ensure_graph_projected(session):
    exists = session.run(
        "CALL gds.graph.exists($name) YIELD exists RETURN exists", name=GRAPH_NAME
    ).single()["exists"]
    if not exists:
        print(f"[gds] projection '{GRAPH_NAME}' absente -> création")
        session.run(
            """
            CALL gds.graph.project(
              $name, 'Protein',
              {INTERACTS_WITH: {orientation: 'UNDIRECTED', properties: 'score'}}
            )
            """,
            name=GRAPH_NAME,
        )


def get_all_toxic_drugs_with_targets(session, min_targets: int) -> pd.DataFrame:
    query = """
    MATCH (d:Drug)-[:TARGETS]->(p:Protein)
    WITH d, count(DISTINCT p) AS n_targets
    WHERE n_targets >= $min_targets
    RETURN d.chembl_id AS chembl_id, d.name AS name, d.warning_class AS warning_class,
           d.hepatox_curated AS hepatox_curated, n_targets
    ORDER BY d.name
    """
    return pd.DataFrame(session.run(query, min_targets=min_targets).data())


def get_seeds(session, chembl_id: str) -> list[str]:
    query = "MATCH (:Drug {chembl_id: $id})-[:TARGETS]->(p:Protein) RETURN DISTINCT p.name AS protein"
    return [r["protein"] for r in session.run(query, id=chembl_id).data()]


def pagerank_top_n(session, seed_names: list[str], top_n: int) -> list[str]:
    query = """
    MATCH (p:Protein) WHERE p.name IN $seed_names
    WITH collect(p) AS seeds
    CALL gds.pageRank.stream($graph, {
      sourceNodes: seeds,
      relationshipWeightProperty: 'score',
      dampingFactor: 0.85
    })
    YIELD nodeId, score
    RETURN gds.util.asNode(nodeId).name AS protein
    ORDER BY score DESC
    LIMIT $top_n
    """
    result = session.run(query, seed_names=seed_names, graph=GRAPH_NAME, top_n=top_n).data()
    return [r["protein"] for r in result]


def induced_subgraph_edges(session, node_names: list[str]) -> nx.Graph:
    query = """
    MATCH (a:Protein)-[r:INTERACTS_WITH]-(b:Protein)
    WHERE a.name IN $nodes AND b.name IN $nodes AND a.name < b.name
    RETURN a.name AS source, b.name AS target, r.score AS score
    """
    G = nx.Graph()
    G.add_nodes_from(node_names)
    for record in session.run(query, nodes=node_names).data():
        G.add_edge(record["source"], record["target"], score=record["score"])
    return G


def motif_nodes_for_drug(session, chembl_id: str, top_n: int, min_clique_size: int) -> set[str]:
    seeds = get_seeds(session, chembl_id)
    top_nodes = pagerank_top_n(session, seeds, top_n)
    if len(top_nodes) < min_clique_size:
        return set()
    G = induced_subgraph_edges(session, top_nodes)
    cliques = [c for c in nx.find_cliques(G) if len(c) >= min_clique_size]
    return set().union(*cliques) if cliques else set()


def jaccard(a: set, b: set) -> float:
    union = a | b
    return len(a & b) / len(union) if union else 0.0


def load_or_compute_motif_nodes(session, drugs: pd.DataFrame, top_n: int, min_clique_size: int) -> dict[str, set]:
    existing_path = COMPARISON_DIR / "drug_motif_nodes.csv"
    cache: dict[str, set] = {}
    if existing_path.exists():
        prev = pd.read_csv(existing_path)
        for _, row in prev.iterrows():
            nodes = row.get("motif_nodes")
            cache[row["chembl_id"]] = set(str(nodes).split("|")) if isinstance(nodes, str) and nodes else set()
        print(f"[cache] {len(cache)} médicaments déjà analysés récupérés depuis {existing_path.name}")

    todo = [cid for cid in drugs["chembl_id"] if cid not in cache]
    print(f"[motifs] {len(todo)} médicaments à calculer (PageRank + cliques, top-{top_n})...")
    new_rows = []
    for chembl_id in tqdm(todo, desc="motifs (nouveaux)"):
        nodes = motif_nodes_for_drug(session, chembl_id, top_n, min_clique_size)
        cache[chembl_id] = nodes
        name = drugs.loc[drugs["chembl_id"] == chembl_id, "name"].iloc[0]
        curated = drugs.loc[drugs["chembl_id"] == chembl_id, "hepatox_curated"].iloc[0]
        new_rows.append({"chembl_id": chembl_id, "name": name, "hepatox_curated": curated,
                          "n_motif_nodes": len(nodes), "motif_nodes": "|".join(sorted(nodes))})

    if new_rows:
        out_path = COMPARISON_DIR / "hepatox_curated_drug_motif_nodes.csv"
        pd.DataFrame(new_rows).to_csv(out_path, index=False)
        print(f"[saved] {out_path} ({len(new_rows)} médicaments nouvellement calculés)")

    return cache


# --------------------------------------------------------------------------
# 5. Test par permutation d'étiquettes (un groupe curé à la fois)
# --------------------------------------------------------------------------
def permutation_test_for_group(
    group_mask: np.ndarray, motif_sets: list[set], rng: np.random.Generator, n_permutations: int
) -> dict:
    n = len(motif_sets)
    idx_a, idx_b = np.triu_indices(n, k=1)
    jac = np.array([jaccard(motif_sets[i], motif_sets[j]) for i, j in zip(idx_a, idx_b)])

    def mean_intra(mask: np.ndarray) -> float:
        ga, gb = mask[idx_a], mask[idx_b]
        same_group = ga & gb
        return jac[same_group].mean() if same_group.any() else np.nan

    observed = mean_intra(group_mask)

    null = np.empty(n_permutations)
    perm_mask = group_mask.copy()
    for p in range(n_permutations):
        rng.shuffle(perm_mask)
        null[p] = mean_intra(perm_mask)

    null_valid = null[~np.isnan(null)]
    p_value = (np.sum(null_valid >= observed) + 1) / (len(null_valid) + 1) if len(null_valid) else np.nan
    z = (observed - null_valid.mean()) / null_valid.std() if len(null_valid) and null_valid.std() > 0 else np.nan

    return {
        "n_drugs_in_group": int(group_mask.sum()),
        "jaccard_intra_group_observed": observed,
        "jaccard_intra_group_null_mean": null_valid.mean() if len(null_valid) else np.nan,
        "jaccard_intra_group_null_std": null_valid.std() if len(null_valid) else np.nan,
        "z_score": z,
        "p_value": p_value,
        "significant_0.05": bool(p_value < 0.05) if not np.isnan(p_value) else False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--list1", type=Path, default=RAW_DIR / "list_Drug_Hepatotoxicity.txt")
    parser.add_argument("--list2", type=Path, default=RAW_DIR / "list_hepa.txt")
    parser.add_argument("--top-n", type=int, default=40, help="Taille du voisinage PageRank (défaut: 40, comme 10_motif_analysis.py)")
    parser.add_argument("--min-clique-size", type=int, default=3)
    parser.add_argument("--min-targets", type=int, default=1)
    parser.add_argument("--pchembl-threshold", type=int, default=PCHEMBL_THRESHOLD)
    parser.add_argument("--n-permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--skip-fetch", action="store_true", help="Ne pas interroger l'API ChEMBL pour les cibles manquantes")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    COMPARISON_DIR.mkdir(parents=True, exist_ok=True)

    if not args.list1.exists() or not args.list2.exists():
        raise FileNotFoundError(
            f"Place les 2 listes curées dans {RAW_DIR} "
            f"(ou passe --list1/--list2), ex: list_Drug_Hepatotoxicity.txt, list_hepa.txt"
        )

    # Avant de lancer, supprime l'ancien fichier de map s'il contient des colonnes erronées
    map_csv = COMPARISON_DIR / "hepatox_curated_chembl_map.csv"

    # --- 1. listes curées + résolution ChEMBL ---
    curated_df = build_curated_sets(args.list1, args.list2)
    resolved_df = resolve_all(curated_df, map_csv)

    unresolved = resolved_df[resolved_df["resolved"] != True][["name", "hepatox_curated"]]
    if len(unresolved):
        unresolved_path = COMPARISON_DIR / "hepatox_curated_unresolved.csv"
        unresolved.to_csv(unresolved_path, index=False)
        print(f"[!] {len(unresolved)} noms non résolus -> {unresolved_path} (revue manuelle, ex. codes de recherche/coquilles)")

    resolved_df = resolved_df[resolved_df["resolved"] == True].copy()

    driver = GraphDatabase.driver(URI, auth=(USER, PASSWORD))
    with driver.session() as session:
        ensure_graph_projected(session)

        # --- 2/3. cibles manquantes + tag Neo4j ---
        if not args.skip_fetch:
            fetch_and_load_missing_targets(session, resolved_df["chembl_id"].tolist(), args.pchembl_threshold)
        tag_curated_drugs(session, resolved_df)

        # --- 4. motifs sur TOUT le fond (curés + reste), avec cache ---
        print(f"\n[setup] médicaments toxiques avec >= {args.min_targets} cible(s)...")
        drugs = get_all_toxic_drugs_with_targets(session, args.min_targets)
        print(f"  -> {len(drugs)} médicaments au total (fond de permutation)")

        motif_cache = load_or_compute_motif_nodes(session, drugs, args.top_n, args.min_clique_size)

    driver.close()

    motif_sets = [motif_cache[cid] for cid in drugs["chembl_id"]]
    hepatox_curated = drugs["hepatox_curated"].to_numpy()

    # --- 5. permutation test : core, puis broad (core U broad-only) ---
    print(f"\n[test] permutation d'étiquettes ({args.n_permutations} tirages), fond = {len(drugs)} médicaments...")
    rows = []
    for group_name, mask_values in [
        ("hepatox_curated_core", hepatox_curated == "core"),
        ("hepatox_curated_broad", np.isin(hepatox_curated, ["core", "broad"])),
    ]:
        res = permutation_test_for_group(mask_values, motif_sets, rng, args.n_permutations)
        res["group"] = group_name
        rows.append(res)

    summary = pd.DataFrame(rows)[
        ["group", "n_drugs_in_group", "jaccard_intra_group_observed", "jaccard_intra_group_null_mean",
         "jaccard_intra_group_null_std", "z_score", "p_value", "significant_0.05"]
    ]
    out_path = COMPARISON_DIR / "hepatox_curated_motif_enrichment.csv"
    summary.to_csv(out_path, index=False)

    print(f"\n{'=' * 88}")
    print("RÉSULTAT -- hépatotoxicité, listes curées (core=intersection, broad=union) vs fond complet")
    print(f"{'=' * 88}")
    print(summary.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print(f"\n[saved] {out_path}")

    # rappel du résultat original (label ChEMBL brut) si disponible
    orig_path = COMPARISON_DIR / "motif_enrichment_by_class.csv"
    if orig_path.exists():
        orig = pd.read_csv(orig_path)
        hep = orig[orig["warning_class"] == "hepatotoxicity"]
        if len(hep):
            print("\n[rappel] résultat original (10_motif_analysis.py, warning_class ChEMBL brut) :")
            print(hep.to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()