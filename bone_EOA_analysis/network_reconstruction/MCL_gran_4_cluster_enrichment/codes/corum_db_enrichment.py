import pandas as pd
from itertools import combinations
from collections import defaultdict
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

# --- Metric Functions ---

def harmonic_mean(a, b):
    return 2 * a * b / (a + b) if (a + b) > 0 else 0

def precision_recall_product(p, b):
    return (len(p & b) ** 2) / (len(p) * len(b)) if p and b else 0

def jaccard(p, b):
    return len(p & b) / len(p | b) if p | b else 0

def semantic_density(group, reference_sets):
    pairs = list(combinations(group, 2))
    if not pairs:
        return 0
    count = sum(1 for x, y in pairs if any(x in s and y in s for s in reference_sets))
    return count / len(pairs)

def compute_weighted_jaccard(P, B):
    jacc_p_vals = [(len(p), max(jaccard(p, b) for b in B)) for p in P]
    jacc_b_vals = [(len(b), max(jaccard(p, b) for p in P)) for b in B]

    jacc_P = sum(size * val for size, val in jacc_p_vals) / sum(size for size, _ in jacc_p_vals)
    jacc_B = sum(size * val for size, val in jacc_b_vals) / sum(size for size, _ in jacc_b_vals)

    return harmonic_mean(jacc_P, jacc_B), jacc_P, jacc_B

# --- Main Pipeline ---

def evaluate_clusters_vs_corum_by_genes(cluster_csv, corum_csv, output_csv):
    # Load cluster file and extract gene names from 'name' column
    clusters_df = pd.read_csv(cluster_csv)
    clusters_df = clusters_df.dropna(subset=["__mcodeCluster", "name"])
    clusters_df["genes"] = clusters_df["name"].astype(str).str.split(";|,")  # Handles semicolon or comma
    clusters_df = clusters_df.explode("genes")
    clusters_df["genes"] = clusters_df["genes"].str.strip()

    clusters = clusters_df.groupby("__mcodeCluster")["genes"].apply(set).to_dict()

    # Load CORUM and parse gene names
    corum_df = pd.read_csv(corum_csv, encoding="ISO-8859-1")
    corum_df = corum_df[corum_df["organism"] == "Human"]
    corum_df["genes"] = corum_df["subunits_gene_name"].fillna("").astype(str).str.split(";").apply(lambda x: set(map(str.strip, x)))
    corum_complexes = corum_df.to_dict("records")

    all_genes = set(clusters_df["genes"].dropna()) | set.union(*(c["genes"] for c in corum_complexes))
    M = len(all_genes)  # total population for hypergeometric

    # Semantic references
    cluster_sets = list(clusters.values())
    corum_sets = [c["genes"] for c in corum_complexes]

    results = []
    pvals = []

    for cluster_id, cluster_genes in clusters.items():
        for corum in corum_complexes:
            complex_genes = corum["genes"]
            overlap = cluster_genes & complex_genes
            if not overlap:
                continue

            # Hypergeometric test (survival function)
            N = len(cluster_genes)  # cluster size
            K = len(complex_genes)  # complex size
            x = len(overlap)        # overlap
            pval = hypergeom.sf(x - 1, M, K, N)

            # Save p-value separately for FDR
            pvals.append(pval)

            # Metrics
            pr = precision_recall_product(cluster_genes, complex_genes)
            jac = jaccard(cluster_genes, complex_genes)
            precision = x / N
            recall = x / K
            fscore = harmonic_mean(precision, recall)
            sem_p = semantic_density(cluster_genes, corum_sets)
            sem_b = semantic_density(complex_genes, cluster_sets)
            semantic = harmonic_mean(sem_p, sem_b)

            results.append({
                "cluster_id": cluster_id,
                "cluster_genes": ";".join(sorted(cluster_genes)),
                "cluster_size": N,
                "complex_size": K,
                "overlap_size": x,
                "overlap_genes": ";".join(sorted(overlap)),
                "F-measure": fscore,
                "Jaccard": jac,
                "PR-product": pr,
                "Semantic Similarity": semantic,
                "Hypergeometric P-value": pval,  # placeholder for now
                "complex_name": corum["complex_name"],
                "organism": corum["organism"],
                "pmid": corum["pmid"],
                "comment_complex": corum["comment_complex"],
                "subunits_uniprot_id": corum["subunits_uniprot_id"],
                "subunits_gene_name": corum["subunits_gene_name"],
                "functions_go_id": corum["functions_go_id"],
            })

    # Apply FDR correction
    _, fdr_corrected, _, _ = multipletests(pvals, method='fdr_bh')
    for i, fdr in enumerate(fdr_corrected):
        results[i]["FDR"] = fdr

    # Save detailed results
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_csv, index=False)
    print(f"✅ Results saved to: {output_csv}")

    # Global Jaccard
    jaccard_total, jacc_P, jacc_B = compute_weighted_jaccard(cluster_sets, corum_sets)
    print(f"📊 Global weighted Jaccard score: {jaccard_total:.4f} (clusters: {jacc_P:.4f}, complexes: {jacc_B:.4f})")

# --- Example usage ---
evaluate_clusters_vs_corum_by_genes(
    r"C:\path\bone_net_MCL_gran_4_all_nodes.csv",  # must-have columns: name (gene), __mcodeCluster
    r"C:\path\corum_humanComplexes_7_2025.csv", # must be downloaded from Corum DB
    r"C:\path\ΗΜ_combs_network_MCODE_DC_2_corumdb_enr.csv"
)
