import pandas as pd
import os

# Define the file path
file_path = r"C:\path_to_DAVID_output.txt"

# Create subfolders for saving DataFrames
output_folder = os.path.join(os.path.dirname(file_path), "clusters")
remade_output_folder = os.path.join(os.path.dirname(file_path), "remade_clusters")
os.makedirs(output_folder, exist_ok=True)
os.makedirs(remade_output_folder, exist_ok=True)

# Read the file line by line
with open(file_path, 'r') as file:
    lines = file.readlines()

# Initialize variables
clusters = {}
scores = {}
current_cluster = None
current_score = None
data = []
columns = None

# Parse the file
for line in lines:
    line = line.strip()
    
    # Identify new annotation clusters
    if line.startswith("Annotation Cluster"):
        if current_cluster is not None and data:
            df = pd.DataFrame(data, columns=columns)
            clusters[f"Cluster_{current_cluster} (Score {current_score})"] = df
            scores[current_cluster] = float(current_score)
            df.to_csv(os.path.join(output_folder, f"Cluster_{current_cluster}_Score_{current_score}.csv"), index=False)
        
        # Extract cluster number and enrichment score
        parts = line.split('\t')
        current_cluster = parts[0].split(" ")[-1]  # Extract cluster number
        current_score = parts[1].split(": ")[-1]  # Extract enrichment score
        data = []  # Reset data list
    elif line.startswith("Category"):
        columns = line.split('\t')  # Extract column names
    else:
        if columns and line:
            data.append(line.split('\t'))

# Save last cluster
if current_cluster is not None and data:
    df = pd.DataFrame(data, columns=columns)
    clusters[f"Cluster_{current_cluster} (Score {current_score})"] = df
    scores[current_cluster] = float(current_score)
    df.to_csv(os.path.join(output_folder, f"Cluster_{current_cluster}_Score_{current_score}.csv"), index=False)

# # Select top 3 clusters with the highest enrichment scores
# top_clusters = sorted(scores, key=scores.get, reverse=True)[:10]

# for cluster in top_clusters:
#     key = f"Cluster_{cluster} (Score {scores[cluster]})"
#     df = clusters[key]
    
#     # Extract relevant data and reformat DataFrame
#     remade_df = pd.DataFrame({
#         "id": range(1, len(df) + 1),
#         "Description": df.iloc[:, 1].str.split("~").str[-1],
#         "proteins": df.iloc[:, 2],  # Extract title without GO ID
#         "q_value": df.iloc[:, 4],  # Set q_value as p-value
#         "p.adjust": df.iloc[:, 4],  # Set p.adjust as p-value
#         "Fold_enrichment": df.iloc[:, -4]  # Assuming the 'Fold Enrichment' column is second-to-last in each cluster's data
#     })
    
#     # Save remade DataFrame
#     remade_df.to_csv(os.path.join(remade_output_folder, f"Remade_Cluster_{cluster}_Score_{scores[cluster]}.csv"), index=False)

# # Print sample output
# for key, df in clusters.items():
#     print(f"\n{key}:")
#     print(df.head())

#keep only mf and bp terms

# Select top clusters with the highest enrichment scores
top_clusters = sorted(scores, key=scores.get, reverse=True)[:10]

for cluster in top_clusters:
    key = f"Cluster_{cluster} (Score {scores[cluster]})"
    df = clusters[key]

    # Filter out CC GO terms, keep only BP and MF
    df_filtered = df[~df['Category'].str.startswith('GOTERM_CC_ALL')].copy()
    
    # Extract relevant data and reformat DataFrame
    remade_df = pd.DataFrame({
        "id": range(1, len(df_filtered) + 1),
        "Description": df_filtered.iloc[:, 1].str.split("~").str[-1],
        "proteins": df_filtered.iloc[:, 2],  # Extract title without GO ID
        "q_value": df_filtered.iloc[:, 4],    # Set q_value as p-value
        "p.adjust": df_filtered.iloc[:, 4],   # Set p.adjust as p-value
        "Fold_enrichment": df_filtered.iloc[:, -4]  # Assuming 'Fold Enrichment' is second-to-last
    })

    # Save remade DataFrame
    remade_df.to_csv(os.path.join(remade_output_folder, f"Remade_Cluster_{cluster}_Score_{scores[cluster]}.csv"), index=False)
