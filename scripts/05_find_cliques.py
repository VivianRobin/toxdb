"""
Extrait l'ego-network d'une protéine cible depuis Neo4j, puis détecte les
cliques maximales avec networkx (Neo4j GDS n'a pas de clique-finder natif).

Usage : python scripts/05_find_cliques.py TP53
"""
import os
import sys

import networkx as nx
from dotenv import load_dotenv
from neo4j import GraphDatabase

load_dotenv()

URI = os.getenv("NEO4J_URI", "bolt://localhost:7687")
USER = os.getenv("NEO4J_USER", "neo4j")
PASSWORD = os.getenv("NEO4J_PASSWORD", "toxkg-dev-pw")

QUERY = """
MATCH (center:Protein {name: $target})-[r:INTERACTS_WITH]-(neighbor)
RETURN center.name AS source, neighbor.name AS target, r.score AS score
"""

# On récupère aussi les arêtes ENTRE voisins (pour que les cliques aient un sens)
QUERY_NEIGHBOR_EDGES = """
MATCH (center:Protein {name: $target})-[:INTERACTS_WITH]-(n1)
MATCH (n1)-[r:INTERACTS_WITH]-(n2)
WHERE (center)-[:INTERACTS_WITH]-(n2)
RETURN n1.name AS source, n2.name AS target, r.score AS score
"""


def build_ego_graph(session, target: str) -> nx.Graph:
    G = nx.Graph()
    G.add_node(target)

    for record in session.run(QUERY, target=target):
        G.add_edge(record["source"], record["target"], score=record["score"])

    for record in session.run(QUERY_NEIGHBOR_EDGES, target=target):
        G.add_edge(record["source"], record["target"], score=record["score"])

    return G


def main():
    if len(sys.argv) != 2:
        print("Usage: python scripts/05_find_cliques.py <NOM_PROTEINE>")
        sys.exit(1)
    target = sys.argv[1]

    driver = GraphDatabase.driver(URI, auth=(USER, PASSWORD))
    with driver.session() as session:
        G = build_ego_graph(session, target)
    driver.close()

    print(f"Ego-network de {target}: {G.number_of_nodes()} noeuds, {G.number_of_edges()} arêtes")

    cliques = list(nx.find_cliques(G))
    cliques.sort(key=len, reverse=True)

    print(f"\n{len(cliques)} cliques maximales trouvées. Top 10 par taille :")
    for c in cliques[:10]:
        print(f"  [{len(c)}] {sorted(c)}")


if __name__ == "__main__":
    main()
