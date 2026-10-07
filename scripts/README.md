# ToxKG — Pipeline Scripts Reference

Numbered, idempotent Python scripts (`scripts/`) building a toxicology
knowledge graph (Neo4j) from STRING (PPI), ChEMBL (drug targets and
withdrawal warnings), and cBioPortal (TCGA PanCancer Atlas mutation and
copy-number data). Each script reads from `data/processed/` (or Neo4j) and
writes new CSVs/figures there, so the pipeline can be re-run from any point
without repeating upstream API calls.

Run order matters: scripts are numbered in dependency order. `01`–`03` build
the graph; `04`+ analyze it. A few later scripts (`13`–`21`) form a second,
independent branch reusing `motif_mutation_crossref.csv` as their common
input rather than chaining strictly 1-by-1 — see notes below.

---

## 01 — `01_download_string.py`
**Input:** none (downloads from `string-db.org`, no API key needed).
**Output:** `data/processed/ppi_edges.csv` (`protein1, protein2,
combined_score`).
Downloads the human (taxid 9606) STRING v12.0 PPI network and filters it to
`combined_score >= STRING_SCORE_THRESHOLD` (default 700, "high confidence").
IDs are still raw STRING IDs at this stage; gene-symbol mapping happens
downstream via the STRING protein-info file.

## 02 — `02_download_chembl_toxic_drugs.py`
**Input:** none (ChEMBL public web service).
**Output:** `data/processed/toxic_drugs.csv` (chembl_id, pref_name,
warning_class, warning_country, warning_year) and
`data/processed/drug_targets.csv` (chembl_id, drug_name, target_chembl_id,
target_gene_symbol, action_type).
Pulls every drug with a ChEMBL `drug_warning` of `warning_type="Withdrawn"`,
plus each drug's official primary target(s) via the `mechanism` endpoint.

## 02b — `02b_download_chembl_offtargets.py`
**Input:** `data/processed/toxic_drugs.csv`.
**Output:** `data/processed/drug_targets_offtarget.csv` (chembl_id,
target_chembl_id, target_gene_symbol, pchembl_value, standard_type).
Adds off-target bioactivity for the same drugs via the ChEMBL `activity`
endpoint, filtered to `pchembl_value >= 6` (sub-µM binding), `assay_type='B'`,
`target_type='SINGLE PROTEIN'`. Captures toxicologically relevant secondary
targets the official mechanism annotation misses (e.g. hERG block by an
antihistamine).

## 02c — `02c_add_drug_by_name.py`
**Input:** drug names as CLI arguments (e.g. `"Terfenadine" "Cisapride"`).
**Output:** merges into `toxic_drugs.csv`, `drug_targets.csv`,
`drug_targets_offtarget.csv` (dedup by chembl_id).
Adds specific drugs by name, independent of whether ChEMBL's `Withdrawn`
flag is populated for them — used to add known cases (e.g. hERG-blocking
drugs) not automatically captured by script 02.

## 03 — `03_load_neo4j.py`
**Input:** `ppi_edges.csv`, `toxic_drugs.csv`, `drug_targets.csv`,
`drug_targets_offtarget.csv` (optional), `gene_mutation_frequency.csv`
(optional, from script 08).
**Output:** populates Neo4j — `(:Protein)`, `(:Protein)-[:INTERACTS_WITH
{score}]-(:Protein)`, `(:Drug {chembl_id, name, warning_class,
warning_year, toxic})`, `(:Drug)-[:TARGETS {action_type, source,
pchembl_value}]->(:Protein)`. Also sets `p.mutation_frequency` when
available.
Loads all processed CSVs into the graph database. Fully idempotent (`MERGE`
everywhere) — safe to re-run after adding new data upstream.

## 04 — `cypher/04_analysis.cypher`
**Input:** the loaded Neo4j graph.
**Output:** none (interactive Cypher queries, not a script).
A set of exploratory Cypher queries for manual inspection of the graph
(not part of the automated pipeline).

## 05 — `05_find_cliques.py`
**Input:** a protein name (CLI arg, e.g. `TP53`); reads from Neo4j.
**Output:** printed maximal cliques (no file).
Extracts a target protein's ego-network (direct neighbors + edges between
them) from Neo4j and finds maximal cliques with `networkx` (Neo4j GDS has no
native clique-finder). Exploratory / single-protein use.

## 06 — `06_compare_drug_networks.py`
**Input:** 2+ drug names (CLI args) and optionally `--reference-genes`;
reads from Neo4j.
**Output:** console summary + `data/processed/comparisons/<drugA>_vs_
<drugB>.csv`.
For each drug: personalized PageRank from its targets (mechanism +
off-target) as seeds, keep top-N. Compares all pairs by Jaccard similarity,
shared/unique proteins, and (optionally) the rank of a reference gene panel.
The manual, pairwise precursor to the systematic analysis in script 10.

## 07 — `07_permutation_test.py`
**Input:** 2 drug names (CLI args); reads from Neo4j.
**Output:** console p-value/z-score + `data/processed/comparisons/
permtest_<A>_vs_<B>.csv` (+ `.png` histogram if matplotlib available).
Tests whether the observed Jaccard between two drugs' PageRank networks is
significant, against a **degree-matched** null (random seed sets drawn from
the same PPI-degree bins as the real seeds — controls for real drug targets
being disproportionately high-degree hub proteins).

## 08 — `08_download_cbioportal_mutations.py`
**Input:** `--from-neo4j` (all proteins currently in the graph) or
`--genes GENE1,GENE2,...`.
**Output:** `data/processed/gene_mutation_frequency.csv` (gene, entrez_id,
n_mutated_samples, n_total_samples, mutation_frequency).
Pulls pan-cancer somatic mutation frequency from the 32 TCGA PanCancer Atlas
studies (whole-exome, discovered dynamically via cBioPortal's `/studies`
endpoint) — chosen over targeted panels (e.g. MSK-IMPACT) to avoid falsely
scoring untested genes as "never mutated."

## 09 — `09_mutation_enrichment_test.py`
**Input:** `--drug NAME --top-n N` (reads from Neo4j) **or** `--genes
GENE1,GENE2,...` (explicit list); requires `gene_mutation_frequency.csv`
loaded into Neo4j (via script 03 after script 08).
**Output:** console p-value/z-score + CSV in `data/processed/comparisons/`.
Tests whether a gene set (a drug's PageRank top-N, or an explicit panel like
the congenital long-QT genes) has significantly higher pan-cancer mutation
frequency than a degree-matched random gene set — same null-model logic as
script 07, applied to mutation frequency instead of network overlap.

## 10 — `10_motif_analysis.py`
**Input:** none beyond the Neo4j graph (CLI flags: `--top-n`,
`--min-clique-size`, `--min-targets`, `--min-drugs-per-class`,
`--n-permutations`).
**Output:** `drug_motif_nodes.csv`, `drug_pair_motif_jaccard.csv`,
`motif_enrichment_by_class.csv` (all in `data/processed/comparisons/`).
The core genome-wide analysis, extending 06/07 from drug pairs to **all**
withdrawn drugs at once. For every drug: PageRank top-N → induced subgraph →
maximal cliques (size ≥ `min-clique-size`) → union = "motif nodes." Computes
Jaccard for all drug pairs in memory, then tests by **label permutation**
(shuffling `warning_class` across drugs, topology fixed) whether drugs
sharing a real toxicity class share more motif nodes than random groupings
of the same size. This produced the paper's headline result: cardiotoxicity
and nephrotoxicity significant after Bonferroni correction; hepatotoxicity
not significant.

## 11 — `11_motif_mutation_crossref.py`
**Input:** `drug_motif_nodes.csv` (script 10), `gene_mutation_frequency.csv`
(script 08), `ppi_edges.csv` (script 01) — all offline, no Neo4j needed.
**Output:** `data/processed/comparisons/motif_mutation_crossref.csv`.
Within each significantly enriched toxicity class, ranks genes by a
composite score combining (a) recurrence across the class's drugs' motifs,
(b) PPI degree (centrality), and (c) pan-cancer mutation frequency. This
file becomes the shared input for most of the downstream scripts (12–21).

## 12 — `12_vizualize_subnetworks.py`
**Input:** `motif_mutation_crossref.csv`, `ppi_edges.csv`.
**Output:** subnetwork figures in `data/processed/visualizations/`.
For each toxicity class, draws the induced PPI subnetwork among its top-N
(`--top-n`, default 15) prioritized genes.

## 13 — `13_enrichment_analysis.py`
**Input:** `motif_mutation_crossref.csv`.
**Output:** enrichment result tables in `data/processed/enrichments/`.
Runs pathway/gene-set enrichment (via `gseapy`, querying Enrichr online:
MSigDB Hallmark 2020 + KEGG 2021 Human) on each class's top-20 prioritized
genes, keeping terms at adjusted p < 0.05.

## 14 — `14_dense_modules.py`
**Input:** `ppi_edges.csv`, `motif_mutation_crossref.csv`.
**Output:** module summary in `data/processed/modules/`.
Builds the PPI graph in `igraph` (edges filtered to score ≥ 400), runs
Louvain/multilevel community detection, and cross-references communities of
size ≥ 5 against the prioritized toxicity genes to find densely connected
modules per class (e.g. the mitochondrial complex-I module for
cardiotoxicity).

## 15 — `15_random_walk.py`
**Input:** `ppi_edges.csv`, `motif_mutation_crossref.csv`.
**Output:** `data/processed/rwr_results/rwr_top_propagated_genes.csv`.
Runs random-walk-with-restart (per toxicity class, seeded from that class's
prioritized genes) on the `igraph`-based PPI network (score ≥ 400) to
produce a propagated top-gene ranking per class — the source of the "top-10
genes per class" used for the CNA burden comparison in the paper.

## 16 — `16_cbioportal_mutations.py`
**Input:** `rwr_results/rwr_top_propagated_genes.csv` (script 15).
**Output:** `data/processed/cbioportal/` (raw per-study CNA data;
`pancancer_detailed_cna.csv` used downstream).
For each class's RWR top-10 genes, fetches GISTIC copy-number alteration
(CNA) values across all TCGA PanCancer Atlas studies via the cBioPortal
REST API directly (`requests`, not the ChEMBL/cBioPortal client library).

## 17 — `17_generate_heatmap.py`
**Input:** `cbioportal/pancancer_detailed_cna.csv` (script 16).
**Output:** `data/processed/cbioportal/toxicity_cna_heatmap.pdf`.
Pivots to a gene × TCGA-study matrix (mean GISTIC value) and renders it as a
`seaborn` heatmap.

## 18 — `18_analyze_top_genes.py`
**Input:** `cbioportal/pancancer_detailed_cna.csv`.
**Output:** `data/processed/cbioportal/top_altered_genes_summary.csv`.
Computes per-gene CNA statistics: fraction of samples with `|GISTIC value|
>= 1` ("altered"), mean absolute GISTIC value, ranked by alteration
frequency. This is the source of the 36.8–37.8% CNA burden figures reported
per class.

## 19 — `19_cross_analysis_toxicity.py`
**Input:** `top_altered_genes_summary.csv` (script 18),
`motif_mutation_crossref.csv` (script 11).
**Output:** `data/processed/cbioportal/genes_with_toxicity_classes.csv` and
`toxicity_class_alteration_summary.csv`.
Merges CNA statistics back onto each gene's toxicity class and aggregates
mean alteration frequency per class — ties the network-motif genes (script
11) to their pan-cancer CNA burden (scripts 16–18) in one table.

## 20 — `20_functionnal_profile.py`
**Input:** `cbioportal/genes_with_toxicity_classes.csv` (script 19).
**Output:** `data/processed/cbioportal/functional_toxicity_profile.csv`.
A short per-class summary (gene count, average alteration %, top gene) —
a compact recap of script 19's output.

## 21 — `21_prioritize_biomarkers.py`
**Input:** `cbioportal/genes_with_toxicity_classes.csv` (script 19).
**Output:** `data/processed/cbioportal/prioritized_biomarkers_panel.csv`.
Computes a class-agnostic composite priority score (`(alteration_freq/100) *
mean_abs_gistic_value`) and keeps the top 10 — the source of the
cross-toxicity candidate biomarker panel (GABARAP, ALOX12, GNB2, GNAS,
GNAI1, PLA2G4A) reported in the paper.

## 22 — `22_hepatox_curated_reanalysis.py`
**Input:** two curated hepatotoxicity drug-name lists (`--list1`, `--list2`,
default `data/raw/hepatox_curated/`); reads/writes Neo4j; optionally queries
the ChEMBL API for drugs not yet in the graph.
**Output:** `data/processed/comparisons/hepatox_curated_chembl_map.csv`,
`hepatox_curated_unresolved.csv`, `hepatox_curated_drug_motif_nodes.csv`,
`hepatox_curated_motif_enrichment.csv`.
Re-runs the script-10 motif/permutation pipeline restricted to a curated,
externally-sourced hepatotoxicity drug set (intersection = "core",
union = "broad" of the two lists) instead of ChEMBL's `warning_class` label,
to test whether hepatotoxicity's null result in script 10 reflects the
noisy/coarse ChEMBL labeling rather than a genuine absence of shared
network motifs. Does not modify `warning_class` or `d.toxic`, so it doesn't
affect scripts 10/11/13/19; adds a separate `d.hepatox_curated` property.
Resilient to ChEMBL API outages: lazy import with retry/backoff, per-name
failure isolation, incremental caching (safe to re-run/resume).

---

## Two analysis branches after script 11

- **12 → 13 → 14 → 15 → 16 → 17 → 18 → 19 → 20 → 21**: the main downstream
  chain used for the published results (subnetwork figures, pathway
  enrichment, community detection, RWR propagation, CNA burden, and the
  final biomarker panel).
- **22**: a standalone, parallel re-analysis of one specific class
  (hepatotoxicity) using better-curated input labels, independent of the
  12–21 chain.

## Common environment

All scripts read Neo4j connection settings from `.env`
(`NEO4J_URI`, `NEO4J_USER`, `NEO4J_PASSWORD`) via `python-dotenv`, and write
into `data/processed/` (raw downloads go to `data/raw/`). See
`requirements.txt` for pinned package versions and `docker-compose.yml` for
the Neo4j (with GDS) service definition.
