import pandas as pd
import networkx as nx
import numpy as np
from scipy.linalg import expm
from node2vec import Node2Vec
from sklearn.metrics import silhouette_score

# --- Load Data ---
# Add here the node info and edge betweenness metrics exported from Cytoscape's Analyze Graph function
# In this example is the Leiden-clustered netwrok. You can input appropriately all other clusterings.
nodes_df = pd.read_csv(r"C:\path\bone_net_Leiden_res_0.1_node_info.csv")
edges_df = pd.read_csv(r"C:\path\network_analytics\bone_net_Leiden_res_0.1_edge_betweenness.csv")

# --- Parse Edge Names ---
edges_df[['source', 'target']] = edges_df['name'].str.extract(r'(\S+)\s+\(interacts with\)\s+(\S+)')

# --- Map Node to Cluster ---
cluster_map = nodes_df.set_index('name')['__leidenCluster'].to_dict()
edges_df['source_cluster'] = edges_df['source'].map(cluster_map)
edges_df['target_cluster'] = edges_df['target'].map(cluster_map)

# --- Label Edges ---
edges_df['edge_type'] = edges_df.apply(
    lambda x: 'intra' if x['source_cluster'] == x['target_cluster'] else 'inter', axis=1
)

# --- Build Graph ---
G = nx.Graph()
for _, row in edges_df.iterrows():
    G.add_edge(row['source'], row['target'])

# Add self-loops if present
for _, row in nodes_df.iterrows():
    if row['SelfLoops'] > 0:
        G.add_edge(row['name'], row['name'])

# ---- GLOBAL NETWORK METRICS ----
global_metrics = {}

# Clustering Coefficient (NetworkX version used only for global metric)
clustering = nx.clustering(G)
global_metrics['Average Clustering Coefficient'] = np.mean(list(clustering.values()))

# Path metrics
if nx.is_connected(G):
    global_metrics['Diameter'] = nx.diameter(G)
    global_metrics['Radius'] = nx.radius(G)
    global_metrics['Average Shortest Path Length'] = nx.average_shortest_path_length(G)
else:
    Gcc = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    global_metrics['Diameter (Largest CC)'] = nx.diameter(Gcc)
    global_metrics['Radius (Largest CC)'] = nx.radius(Gcc)
    global_metrics['Average Shortest Path Length (Largest CC)'] = nx.average_shortest_path_length(Gcc)

# # Centralities (retain only unique ones)
# eig_cent = nx.eigenvector_centrality_numpy(G)
# A = nx.to_numpy_array(G)
# subgraph_cent = np.diag(expm(A))

# ---- INTRA-CLUSTER METRICS (using Cytoscape node attributes) ----
cyto_metrics = ['Degree', 'ClusteringCoefficient', 'BetweennessCentrality', 'ClosenessCentrality', 'Eccentricity']
intra_cluster_metrics = {}

for cluster in nodes_df['__leidenCluster'].dropna().unique():
    cluster_nodes = nodes_df[nodes_df['__leidenCluster'] == cluster]
    node_names = cluster_nodes['name'].tolist()
    subgraph = G.subgraph(node_names)

    if subgraph.number_of_nodes() == 0:
        continue  # Skip empty clusters

    metrics = {
        'Num Nodes': subgraph.number_of_nodes(),
        'Num Edges': subgraph.number_of_edges(),
        'Density': nx.density(subgraph),
    }

    # Path length if connected
    if nx.is_connected(subgraph):
        metrics['Avg Shortest Path Length'] = nx.average_shortest_path_length(subgraph)

    for metric in cyto_metrics:
        if metric in cluster_nodes.columns:
            metrics[f"Avg {metric} (Cytoscape)"] = cluster_nodes[metric].mean()
            metrics[f"Std {metric} (Cytoscape)"] = cluster_nodes[metric].std()

    metrics['Num Isolated Nodes'] = cluster_nodes['IsSingleNode'].sum()
    intra_cluster_metrics[int(cluster)] = metrics

# --- AGGREGATED CLUSTER METRICS ---
cluster_density_vals = [m['Density'] for m in intra_cluster_metrics.values()]
cluster_degree_vals = [m['Avg Degree (Cytoscape)'] for m in intra_cluster_metrics.values() if 'Avg Degree (Cytoscape)' in m]
cluster_clust_vals = [m['Avg ClusteringCoefficient (Cytoscape)'] for m in intra_cluster_metrics.values() if 'Avg ClusteringCoefficient (Cytoscape)' in m]
cluster_sizes = [m['Num Nodes'] for m in intra_cluster_metrics.values()]

# ---- INTER-CLUSTER METRICS ----
inter_edges = edges_df[edges_df['edge_type'] == 'inter']
inter_cluster_matrix = pd.crosstab(inter_edges['source_cluster'], inter_edges['target_cluster'])

bridge_nodes = [n for n in G.nodes() if len(set(
    cluster_map.get(nb) for nb in G[n] if cluster_map.get(nb) != cluster_map.get(n)
)) > 0]

# --- NODE2VEC EMBEDDING & SILHOUETTE (Using Cytoscape cluster labels) ---

# Keep only nodes with a known cluster
valid_nodes = [n for n in G.nodes() if not pd.isna(cluster_map.get(n))]
node_labels = [cluster_map[n] for n in valid_nodes]

# Generate Node2Vec embeddings
node2vec = Node2Vec(G.subgraph(valid_nodes), dimensions=64, walk_length=30, num_walks=200, workers=2, seed=42)
model = node2vec.fit(window=10, min_count=1)

# Get embeddings in correct node order
embeddings = [model.wv[str(n)] for n in valid_nodes]

# Compute silhouette score using Cytoscape's cluster labels
if len(set(node_labels)) > 1:
    silhouette_emb_score = silhouette_score(embeddings, node_labels)
else:
    silhouette_emb_score = float("nan")

# Store or print
global_metrics['Silhouette Coefficient (Node2Vec + Cytoscape clusters)'] = silhouette_emb_score

# ---- EXPORT RESULTS ----
with open(r"C:\path\network_analytics\Leiden_res_0.1_metrics.txt", "w") as f:
    # f.write("=== GLOBAL NETWORK METRICS ===\n")
    # for k, v in global_metrics.items():
    #     f.write(f"{k}: {v}\n")

    # f.write("\n=== CENTRALITY METRICS (Averages) ===\n")
    # f.write(f"Avg Eigenvector Centrality: {np.mean(list(eig_cent.values()))}\n")
    # f.write(f"Avg Subgraph Centrality: {np.mean(subgraph_cent)}\n")

    f.write("\n=== INTRA-CLUSTER METRICS ===\n")
    for c, m in intra_cluster_metrics.items():
        f.write(f"\nCluster {c}:\n")
        for k, v in m.items():
            f.write(f"  {k}: {v}\n")

    # f.write("\n=== INTER-CLUSTER METRICS ===\n")
    # f.write(f"Number of Inter-cluster Edges: {len(inter_edges)}\n")
    # f.write(f"Number of Bridge Nodes: {len(bridge_nodes)}\n")
    # f.write("\nInter-Cluster Edge Matrix:\n")
    # f.write(inter_cluster_matrix.to_string())

    f.write("\n=== AGGREGATED CLUSTER SUMMARY ===\n")
    f.write(f"Avg Cluster Density: {np.mean(cluster_density_vals)}\n")
    f.write(f"Avg Cluster Clustering Coefficient (Cytoscape): {np.mean(cluster_clust_vals)}\n")
    f.write(f"Avg Cluster Degree (Cytoscape): {np.mean(cluster_degree_vals)}\n")
    f.write(f"Max Cluster Size: {max(cluster_sizes)}\n")
    f.write(f"Min Cluster Size: {min(cluster_sizes)}\n")
    f.write(f"Cluster Size Std Dev: {np.std(cluster_sizes)}\n")
    f.write(f"Cluster Degree Heterogeneity (Cytoscape): {np.std(cluster_degree_vals)}\n")
    f.write(f"Modularity Proxy Score: {np.mean(cluster_density_vals) * len(cluster_sizes)}\n")

    f.write("\n=== Silhouette Coefficient (Node2Vec + Cytoscape Clusters) ===\n")
    f.write(f"Silhouette Coefficient: {global_metrics['Silhouette Coefficient (Node2Vec + Cytoscape clusters)']}\n")

