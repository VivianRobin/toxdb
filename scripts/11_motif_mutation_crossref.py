"""
Croise `drug_motif_nodes.csv` (script 10) avec `gene_mutation_frequency.csv`
(script 08) pour prioriser, au sein des classes de toxicite significativement
enrichies en motifs partages (script 10 -> motif_enrichment_by_class.csv),
les genes qui sont a la fois :

  (a) topologiquement recurrents : presents dans le "motif" (clique commune)
      d'une large fraction des medicaments de la classe
  (b) topologiquement centraux dans le reseau IPP global : degre eleve
      (nombre de partenaires directs, calcule depuis ppi_edges.csv)
  (c) frequemment mutes en pan-cancer (gene_mutation_frequency.csv, TCGA)

Ne necessite PAS Neo4j : tout est recalcule depuis les CSV deja produits
par les scripts 08 et 10, donc executable hors ligne.

Usage :
  python scripts/11_motif_mutation_crossref.py
  python scripts/11_motif_mutation_crossref.py --alpha 0.05 --top-n 15
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

PROC_DIR = Path(__file__).resolve().parents[1] / "data" / "processed"
COMPARISON_DIR = PROC_DIR / "comparisons"


def load_ppi_degree(ppi_path: Path) -> pd.Series:
    """Degre non-oriente de chaque proteine dans le reseau IPP (STRING)."""
    edges = pd.read_csv(ppi_path)
    deg = pd.concat([edges["source"], edges["target"]]).value_counts()
    deg.name = "ppi_degree"
    return deg


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--alpha", type=float, default=0.05, help="Seuil de significativite (defaut: 0.05)")
    parser.add_argument("--top-n", type=int, default=15, help="Nb de genes prioritaires affiches par classe")
    parser.add_argument(
        "--min-recurrence",
        type=float,
        default=0.0,
        help="Filtre optionnel : fraction min. de medicaments de la classe portant le gene (0-1)",
    )
    args = parser.parse_args()

    enrich = pd.read_csv(COMPARISON_DIR / "motif_enrichment_by_class.csv")
    sig_col = [c for c in enrich.columns if c.startswith("significant_")][0]
    sig_classes = enrich.loc[enrich[sig_col], "warning_class"].tolist()
    print(f"[classes significatives, alpha={args.alpha}] {sig_classes}")

    motif_nodes = pd.read_csv(COMPARISON_DIR / "drug_motif_nodes.csv")
    mutfreq = pd.read_csv(PROC_DIR / "gene_mutation_frequency.csv")
    mutfreq_lookup = mutfreq.set_index("gene")[["mutation_frequency", "n_mutated_samples", "n_total_samples"]]

    ppi_deg = load_ppi_degree(PROC_DIR / "ppi_edges.csv")

    all_results = []

    for cls in sig_classes:
        sub = motif_nodes[motif_nodes["warning_class"] == cls]
        n_drugs_cls = sub["chembl_id"].nunique()
        if n_drugs_cls == 0:
            continue

        # frequence de recurrence de chaque gene dans les motifs de la classe
        gene_drug_count = {}
        for nodes_str in sub["motif_nodes"].dropna():
            for g in set(nodes_str.split("|")):
                gene_drug_count[g] = gene_drug_count.get(g, 0) + 1

        df = pd.DataFrame(
            {"gene": list(gene_drug_count.keys()), "n_drugs_with_motif": list(gene_drug_count.values())}
        )
        df["n_drugs_in_class"] = n_drugs_cls
        df["recurrence_frac"] = df["n_drugs_with_motif"] / n_drugs_cls
        df["warning_class"] = cls

        df = df.merge(mutfreq_lookup, left_on="gene", right_index=True, how="left")
        df = df.merge(ppi_deg, left_on="gene", right_index=True, how="left")
        df["mutation_frequency"] = df["mutation_frequency"].fillna(0.0)
        df["ppi_degree"] = df["ppi_degree"].fillna(0).astype(int)

        if args.min_recurrence > 0:
            df = df[df["recurrence_frac"] >= args.min_recurrence]

        # --------------------------------------------------------------------------------
        # NOUVELLE MÉTHODE : Correction du biais topologique par Modèle Nul (Binning)
        # --------------------------------------------------------------------------------
        
        # 1. Sécurité : adapter le nombre de bins au nombre de gènes disponibles
        n_bins = min(10, len(df))
        
        if n_bins > 1:
            # Création des groupes (bins) basés sur le degré IPP
            df['degree_bin'] = pd.qcut(df['ppi_degree'].rank(method='first'), q=n_bins, labels=False)
            
            # Calcul de la récurrence attendue (R_null) pour chaque groupe de degré
            null_expectation = df.groupby('degree_bin')['recurrence_frac'].mean().reset_index()
            null_expectation.rename(columns={'recurrence_frac': 'R_null'}, inplace=True)
            
            # Fusion et calcul de l'enrichissement par rapport au bruit de fond
            df = df.merge(null_expectation, on='degree_bin', how='left')
            epsilon = 0.05 
            df['enrichment_score'] = (df['recurrence_frac'] - df['R_null']) / (df['R_null'] + epsilon)
        else:
            # Cas limite : s'il y a 1 seul gène, l'enrichissement est juste sa récurrence
            df['enrichment_score'] = df['recurrence_frac']
            
        # 2. Intégration de la mutation comme un bonus pondéré (+20% max)
        df["rank_mutation_freq"] = df["mutation_frequency"].rank(pct=True)
        
        # Un score d'enrichissement négatif signifie que le gène est moins récurrent qu'attendu :
        # on clip à 0 pour ne garder que les gènes vraiment enrichis avant d'appliquer le bonus.
        df['enrichment_score'] = df['enrichment_score'].clip(lower=0) 
        df['priority_score'] = df['enrichment_score'] * (1 + (df['rank_mutation_freq'] * 0.2))

        # 3. Nettoyage des colonnes de calcul temporaires
        cols_to_drop = ['degree_bin', 'R_null', 'enrichment_score']
        df.drop(columns=[c for c in cols_to_drop if c in df.columns], inplace=True)
        # --------------------------------------------------------------------------------

        df = df.sort_values("priority_score", ascending=False)
        all_results.append(df)

        print(f"\n{'=' * 90}")
        print(f"TOP {args.top_n} genes prioritaires -- {cls}  ({n_drugs_cls} medicaments)")
        print(f"{'=' * 90}")
        cols = ["gene", "recurrence_frac", "ppi_degree", "mutation_frequency", "n_mutated_samples", "priority_score"]
        print(df[cols].head(args.top_n).to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    if not all_results:
        print("[warning] aucune classe significative -- rien a croiser")
        return

    combined = pd.concat(all_results, ignore_index=True)
    out_cols = [
        "warning_class", "gene", "n_drugs_with_motif", "n_drugs_in_class", "recurrence_frac",
        "ppi_degree", "mutation_frequency", "n_mutated_samples", "n_total_samples", "priority_score",
    ]
    out_path = COMPARISON_DIR / "motif_mutation_crossref.csv"
    combined[out_cols].to_csv(out_path, index=False)
    print(f"\n[saved] {out_path}  ({len(combined)} lignes, {combined['warning_class'].nunique()} classes)")


if __name__ == "__main__":
    main()
