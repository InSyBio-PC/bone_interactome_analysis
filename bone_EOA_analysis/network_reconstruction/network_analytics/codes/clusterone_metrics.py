import pandas as pd
import networkx as nx
from sklearn.metrics import silhouette_score
import numpy as np

# --- Step 1: Load the cluster file ---
cluster_df = pd.read_csv(r"C:\path\network_analytics\bone_ppi_clusterone_results_size_4_dens_0.3.csv")  # Adjust separator if needed. Exported from ClusterONE plugin in Cytoscape
cluster_df['Members'] = cluster_df['Members'].apply(lambda x: x.split())

# Build gene -> cluster label map
gene_to_cluster = {}
for idx, row in cluster_df.iterrows():
    for gene in row['Members']:
        gene_to_cluster[gene] = row['Cluster']

# --- Step 2: Load the edge file ---
edge_df = pd.read_csv(r"C:\path\network_analytics\final_bone_ppi_network_file_default_edge.csv")
edge_pairs = edge_df['name'].str.extract(r'(\w+)\s+\(interacts with\)\s+(\w+)')
edges = list(edge_pairs.itertuples(index=False, name=None))

# --- Step 3: Build the graph ---
G = nx.Graph()
G.add_edges_from(edges)

# Make sure all cluster genes are in the graph (add isolated nodes if needed)
all_cluster_genes = set(gene_to_cluster.keys())
G.add_nodes_from(all_cluster_genes)

# --- Step 4: Compute metrics ---

# Clustering Coefficient
clustering_coeff = nx.average_clustering(G)

# Density
density = nx.density(G)

# Heterogeneity: Coefficient of Variation of node degrees
degrees = [deg for _, deg in G.degree()]
heterogeneity = np.std(degrees) / np.mean(degrees) if np.mean(degrees) != 0 else 0

# Average Degree
avg_degree = np.mean(degrees)

# Silhouette Coefficient: Requires distance matrix + cluster labels
nodes_in_clusters = list(gene_to_cluster.keys())
node_idx = {node: i for i, node in enumerate(nodes_in_clusters)}

# Build distance matrix from shortest paths
lengths = dict(nx.all_pairs_shortest_path_length(G))
n = len(nodes_in_clusters)
dist_matrix = np.full((n, n), np.inf)

for src, targets in lengths.items():
    if src in node_idx:
        for tgt, dist in targets.items():
            if tgt in node_idx:
                dist_matrix[node_idx[src], node_idx[tgt]] = dist

# Replace infs (disconnected pairs) with a large constant
max_dist = np.nanmax(dist_matrix[dist_matrix != np.inf])
dist_matrix[dist_matrix == np.inf] = max_dist + 1

# Cluster labels for silhouette score
labels = [gene_to_cluster[node] for node in nodes_in_clusters]

# Compute silhouette score (using precomputed distance matrix)
silhouette = silhouette_score(dist_matrix, labels, metric='precomputed')

# --- Step 5: Print results ---
print("=== Network-Based Clustering Metrics ===")
print(f"Average Clustering Coefficient: {clustering_coeff:.4f}")
print(f"Graph Density: {density:.4f}")
print(f"Heterogeneity (CV of Degree): {heterogeneity:.4f}")
print(f"Average Node Degree: {avg_degree:.2f}")
print(f"Silhouette Coefficient: {silhouette:.4f}")
