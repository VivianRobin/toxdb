// ============================================================
// Analyses reproduisant la méthodo de la thèse, en Cypher + GDS
// A exécuter dans Neo4j Browser (http://localhost:7474) ou via
// le driver Python, une fois 03_load_neo4j.py passé.
// ============================================================

// --- 0. Vue d'ensemble ---
MATCH (p:Protein) RETURN count(p) AS n_proteins;
MATCH ()-[r:INTERACTS_WITH]-() RETURN count(r) AS n_interactions;
MATCH (d:Drug) RETURN count(d) AS n_toxic_drugs;
MATCH (:Drug)-[t:TARGETS]->(:Protein) RETURN count(t) AS n_drug_target_links;


// --- 1. Projeter le graphe dans GDS (nécessaire pour tous les algos) ---
CALL gds.graph.project(
  'ppi',
  'Protein',
  {INTERACTS_WITH: {orientation: 'UNDIRECTED', properties: 'score'}}
);


// --- 2. Ego-network d'une protéine cible (équivalent Objectif 1, étape "ego network") ---
// Remplace 'TP53' par la protéine cible du principe actif étudié
MATCH (center:Protein {name: 'TP53'})-[:INTERACTS_WITH]-(neighbor)
RETURN center.name AS center, collect(neighbor.name) AS ego_network_1hop;

// Version 2 sauts (voisins des voisins) si le 1-hop est trop petit
MATCH (center:Protein {name: 'TP53'})-[:INTERACTS_WITH*1..2]-(neighbor)
RETURN center.name AS center, collect(DISTINCT neighbor.name) AS ego_network_2hop;


// --- 3. RWR (Random Walk with Restart) depuis les cibles d'un médicament toxique ---
// équivalent de "RWR et filtration du classement des n protéines associées
// à une protéine cible du principe actif" (étapes clés Objectif 1)
MATCH (d:Drug {name: 'NOM_DU_MEDICAMENT'})-[:TARGETS]->(seed:Protein)
WITH collect(seed) AS seeds
CALL gds.pageRank.stream('ppi', {
  sourceNodes: seeds,
  relationshipWeightProperty: 'score',
  dampingFactor: 0.85
})
YIELD nodeId, score
RETURN gds.util.asNode(nodeId).name AS protein, score
ORDER BY score DESC
LIMIT 50;


// --- 4. Détection de communautés (= "arrangements communs d'interactions") ---
// Louvain sur l'ego-network élargi, pour faire émerger des sous-modules
// équivalents à tes "arrangements communs IPP"
CALL gds.louvain.stream('ppi', {relationshipWeightProperty: 'score'})
YIELD nodeId, communityId
RETURN communityId, collect(gds.util.asNode(nodeId).name) AS members
ORDER BY size(members) DESC
LIMIT 20;


// --- 5. Cliques (sous-groupes totalement connectés) dans un ego-network ---
// GDS n'a pas de clique-finder natif -> on exporte le sous-graphe et on
// utilise networkx côté Python (voir scripts/05_find_cliques.py)
MATCH (center:Protein {name: 'TP53'})-[r:INTERACTS_WITH]-(neighbor)
RETURN center.name AS source, neighbor.name AS target, r.score AS score;


// --- 6. Comparaison de deux ego-networks (= "alignement de graphe" simplifié) ---
// Similarité de Jaccard entre les voisinages de deux protéines cibles de
// deux principes actifs différents -- première approximation avant un vrai
// alignement de graphe
MATCH (a:Protein {name: 'PROTEINE_CIBLE_1'})-[:INTERACTS_WITH]-(na)
WITH a, collect(na.name) AS neighborsA
MATCH (b:Protein {name: 'PROTEINE_CIBLE_2'})-[:INTERACTS_WITH]-(nb)
WITH a, neighborsA, b, collect(nb.name) AS neighborsB
WITH a, b, neighborsA, neighborsB,
     [x IN neighborsA WHERE x IN neighborsB] AS intersection
RETURN a.name, b.name,
       size(intersection) AS shared_neighbors,
       toFloat(size(intersection)) / size(apoc.coll.toSet(neighborsA + neighborsB)) AS jaccard_similarity;


// --- Nettoyage de la projection GDS quand t'as fini ---
CALL gds.graph.drop('ppi');
