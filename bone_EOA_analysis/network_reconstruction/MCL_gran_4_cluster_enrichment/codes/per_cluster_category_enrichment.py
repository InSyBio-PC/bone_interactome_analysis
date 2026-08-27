import pandas as pd
import numpy as np
import requests
from scipy.stats import hypergeom
import os

# === Load Files ===
cluster_file = pd.read_csv(r"C:\path\bone_net_MCL_gran_4_all_nodes.csv")  # must-have columns: name (gene), __mcodeCluster
main_cat_map = pd.read_excel(r"C:\path\Protein_Groups.xlsx")  # Gene name, Main Category
uniprot_map = pd.read_csv(
    r"C:\path\UNIPROT_PPI_MAPPING\uniprotkb_HUMAN_REVIEWED_GENES_REFSEQ_SEQUENCES_2025_07_01.tsv", # Downloaded from offical UniProt's resources
    sep="\t"
)  # Entry, Gene Names

# === Build gene -> UID map (first match only) ===
gene_to_uid = {}
for _, row in uniprot_map.iterrows():
    uid = row["Entry"]
    gene_names = str(row["Gene Names"]).split()
    for gene in gene_names:
        if gene not in gene_to_uid:
            gene_to_uid[gene] = uid

# === Build gene -> Main Category map from provided mapping ===
gene_to_cat = dict(zip(main_cat_map["Gene name"], main_cat_map["Main Category"]))

# === Map clustering file genes to UID and Main Category ===
cluster_file.rename(columns={"name": "gene"}, inplace=True)
cluster_file["uid"] = cluster_file["gene"].map(gene_to_uid)
cluster_file["Main Category"] = cluster_file["gene"].map(gene_to_cat)
cluster_file["details"] = None
cluster_file["API_check"] = 0

# === Define inference logic (unchanged) ===
def infer_protein_group(data):
    details = []
    cell_loc_evidence = []

    try:
        name = data["proteinDescription"]["recommendedName"]["fullName"]["value"].lower()
    except:
        name = ""

    keywords = [kw.get("value", "").lower() for kw in data.get("keywords", [])]
    locations = []
    for comment in data.get("comments", []):
        if comment.get("commentType") == "SUBCELLULAR LOCATION":
            for loc in comment.get("subcellularLocations", []):
                loc_str = loc.get("location", {}).get("value", "")
                locations.append(loc_str.lower())
                cell_loc_evidence.append(loc_str)

    if "apolipoprotein" in name:
        return "Apolipoproteins", "; ".join(["protein name: apolipoprotein"] + cell_loc_evidence)
    if "immunoglobulin" in keywords or "immunoglobulin" in name:
        return "Immunoglobulins", "; ".join(["keyword or name: immunoglobulin"] + cell_loc_evidence)
    if "milk" in name:
        return "Milk-specific proteins", "; ".join(["protein name: milk"] + cell_loc_evidence)
    if "membrane" in locations:
        return "Membrane_Proteins", "; ".join(["subcellularLocation: membrane"] + cell_loc_evidence)
    if "extracellular matrix" in locations or "extracellular" in locations:
        if "structural" in name:
            return "Structural_ECM_proteins", "; ".join(["name: structural ECM"] + cell_loc_evidence)
        elif "enzyme" in name or "peptidase" in name:
            return "ECM_Enzymes", "; ".join(["name: ECM enzyme"] + cell_loc_evidence)
        elif "matricellular" in name:
            return "Matricellular", "; ".join(["name: matricellular"] + cell_loc_evidence)
        else:
            return "Other ECM proteins", "; ".join(["location: ECM (other)"] + cell_loc_evidence)
    if "secreted" in keywords or "secreted" in locations:
        return "Secreted_Factors", "; ".join(["keyword or location: secreted"] + cell_loc_evidence)
    if "serum" in name:
        return "Serum_Proteins", "; ".join(["name: serum"] + cell_loc_evidence)
    if "enzyme" in name or "kinase" in name or "phosphatase" in name:
        return "Cellular Enzymes", "; ".join(["name: enzyme"] + cell_loc_evidence)
    if name:
        return "Other_Cellular_Proteins", "; ".join(["fallback: default to other cellular"] + cell_loc_evidence)
    return "Other_Cellular_Proteins", "; ".join(cell_loc_evidence)

# === Infer missing Main Categories using API ===
for idx, row in cluster_file[cluster_file["Main Category"].isna()].iterrows():
    uid = row["uid"]
    if pd.isna(uid):
        cluster_file.at[idx, "details"] = "No UID available"
        continue
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.json"
    try:
        r = requests.get(url)
        if r.status_code == 200:
            data = r.json()
            category, detail = infer_protein_group(data)
            cluster_file.at[idx, "Main Category"] = category
            cluster_file.at[idx, "details"] = detail
            cluster_file.at[idx, "API_check"] = 1
        else:
            cluster_file.at[idx, "details"] = f"Failed to fetch: {r.status_code}"
    except Exception as e:
        cluster_file.at[idx, "details"] = f"Error: {str(e)}"

# === Functional Enrichment Analysis per cluster ===
output_dir = r"C:\Users\harry\Desktop\protein_group__enrichment_results"
os.makedirs(output_dir, exist_ok=True)

# Background (population): all genes in cluster_file with category info
background = cluster_file.dropna(subset=["Main Category"])
M = background["gene"].nunique()  # total population size
background_cat_counts = background["Main Category"].value_counts().to_dict()

for cluster_id, subdf in cluster_file.dropna(subset=["__mcodeCluster"]).groupby("__mcodeCluster"):
    subdf = subdf.dropna(subset=["Main Category"])
    genes_in_cluster = subdf["gene"].unique()
    N = len(genes_in_cluster)
    cluster_cat_counts = subdf["Main Category"].value_counts()

    enrichment_results = []
    for cat, x in cluster_cat_counts.items():
        K = background_cat_counts.get(cat, 0)
        # Hypergeometric test: P(X>=x)
        pval = hypergeom.sf(x - 1, M, K, N)
        fold_enrichment = (x / N) / (K / M) if K > 0 else 0
        enrichment_results.append({
            "Term": cat,
            "description": cat,
            "Category": cat,
            "Count": x,
            "%": round(x / N * 100, 2),
            "PValue": pval,
            "Genes": "; ".join(subdf[subdf["Main Category"] == cat]["gene"]),
            "List Total": N,
            "Pop Total": M,
            "Fold Enrichment": round(fold_enrichment, 2),
        })

    # Adjust P-values using Benjamini-Hochberg (FDR)
    df_enrich = pd.DataFrame(enrichment_results).sort_values("PValue").reset_index(drop=True)
    m = len(df_enrich)
    df_enrich["FDR"] = df_enrich["PValue"] * m / (df_enrich.index + 1)
    df_enrich["FDR"] = df_enrich["FDR"].clip(upper=1.0)

    # Save results
    cluster_file_safe = str(cluster_id).replace(" ", "_")
    df_enrich.to_csv(os.path.join(output_dir, f"cluster_{cluster_file_safe}_enrichment.csv"), index=False)

print("✅ Enrichment analysis complete. Results saved to:", output_dir)
