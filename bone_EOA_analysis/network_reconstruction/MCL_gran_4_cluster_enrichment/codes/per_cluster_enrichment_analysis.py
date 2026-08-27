#----------------------------------------------------------------------------------------------------------------#
# Cytoscape per Cluster GO Enrichment Analysis
#----------------------------------------------------------------------------------------------------------------#

"""
GO Enrichment Analysis and Circular Plotting per Cluster

This script performs Gene Ontology (GO) enrichment analysis on gene clusters using
gProfiler with two background options (human genome and custom cluster genes).

!!IMPORTANT!! As a default, custom cluster genes is the background and contains the total genes from the Cytoscape newtork.
If you want to find the enrichment of your clusters against a custom list, update the 'background_mode' variable.
Also you must provide a csv with a single column named 'gene' containing the desired background genes. 
SOS! CORUM DB and Protein Groups Enrichment use ONLY the set of network's genes as a background

The script generates circular chord plots visualizing enriched GO terms and associated genes,
plus a neat legend PNG explaining GO term colors and categories.

User provides:
- Input CSV path with Cytoscape luster info:
File must be the exported Cytoscape Clustered network node table. 
!!IMPORTANT!! For different clustering algorithms set the appropriate value on the 'clustering_type' variable
Values:
> __mcodeCluster (If MCODE was used)
> __mclCluster (If MCL was used) 
> __leidenCluster (If Leiden was used) 
> etc. based on your clustering node table. Generally it is __XCluster, where X is the lowercased algorithm name

! You can input your own csv with clustering results as long as it contains two columns: 
1) the genes of the network in a 'name' column, and 
2) the clusters they are assigned to (e.g. gene1 --> 2, gene2 --> 5, etc) in a column named X. !!But you have to set X in the 'clustering_type' variable!!
! You can also perform gene enrichment analysis on an unclustered list of genes against the human genome background or a list of your choice
Just provide a csv and set a __XCluster column of your choice with only one value (e.g. 1) and all your genes in a 'name' column.
The rest follows as normally.

- Output base directory for results and plots

Outputs:
- Enrichment CSV files per cluster and background (human and custom which can be either the network genes or a custom background)
- Circular chord plots per cluster for the enriched GO terms
- GO term legends as PNG files per cluster
- Protein Group & CORUM DB complex Enrichment results and plots per cluster (if group_enrichment_mode=1, set to 0 if you don't want it)
SOS! 
- For Protein Group Mapping you also need to download "uniprotkb_HUMAN_REVIEWED_GENES_REFSEQ_SEQUENCES_2025_07_01.tsv" from UniProt DB's resources for appropriate mapping
- CORUM DB and Protein Groups Enrichment use ONLY the set of network's genes as a background


"""

import os
import ast
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from gprofiler import GProfiler
import requests
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
from pathlib import Path

#----------------------------------------------------------------------------------------------------------------#
# Set directories
#----------------------------------------------------------------------------------------------------------------#
#DATASET_DIR=r"C:\path\Enrichment_Analysis
DATASET_DIR=Path(__file__).resolve().parent
print(DATASET_DIR)
#Set input csv path containing genes and clusters
CSV_PATH = r"C:\path\bone_net_MCL_gran_4_all_nodes.csv"
#Set output csv path
OUTPUT_BASE = r"C:\path\MCL_gran_4_cluster_enrichment\total"

BACKGROUND_DIRS = [
    os.path.join(OUTPUT_BASE, "human_background"),
    os.path.join(OUTPUT_BASE, "custom_background")
]
#----------------------------------------------------------------------------------------------------------------#
# Set clusterng type:
# > "__mcodeCluster" (If MCODE was used)
# > "__mclCluster" (If MCL was used) 
# > "__leidenCluster" (If Leiden was used) 
# > etc. based on your clustering node table.
#----------------------------------------------------------------------------------------------------------------#
clustering_type= '__mclCluster'
#----------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------------#
# Set enrichment background mode:
# > "default"     (uses the set of genes from the loaded network as a background)
# > "custom_list" (use external gene list CSV with a 'gene' column)
#----------------------------------------------------------------------------------------------------------------#
background_mode = "default"
EXTERNAL_GENE_LIST_PATH = r"C:\path\background.csv"  # Path to CSV with 'gene' column


#----------------------------------------------------------------------------------------------------------------#
# Perform Protein Group & CORUM DB complex Enrichment per cluster
# > Set to 1 if you want, 0 if you don't 
#----------------------------------------------------------------------------------------------------------------#
group_enrichment_mode = 1
#----------------------------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------------------------#
# PERFORM GO ENRICHMENT PER CLUSTER
#----------------------------------------------------------------------------------------------------------------#

# Map gProfiler 'source' abbreviations to full Category strings
SOURCE_CATEGORY_MAP = {
    'GO:BP': 'GOTERM_BP_DIRECT', 'GO:MF': 'GOTERM_MF_DIRECT', 'GO:CC': 'GOTERM_CC_DIRECT',
    'KEGG': 'KEGG_PATHWAY', 'REAC': 'REACTOME_PATHWAY', 'WP': 'WIKIPATHWAYS',
    'TF': 'TRANSFAC', 'MIRNA': 'MIRNA_TARGETS', 'HPA': 'HUMAN_PROTEIN_ATLAS',
    'CORUM': 'CORUM_COMPLEX', 'HP': 'HUMAN_PHENOTYPE_ONTOLOGY'
}
GO_CATEGORIES = ['GOTERM_BP_DIRECT', 'GOTERM_CC_DIRECT', 'GOTERM_MF_DIRECT']

# === Functions ===
def load_gene_df(path):
    df = pd.read_csv(path)
    return df.dropna(subset=['name', clustering_type])

def parse_genes(x):
    if isinstance(x, str) and x.startswith('['):
        try: val = ast.literal_eval(x)
        except: val = []
    else:
        val = x
    if isinstance(val, list):
        flat = []
        for i in val:
            flat.extend(i if isinstance(i, list) else [i])
        return ','.join(flat)
    return ''

def load_external_gene_list(path):
    df = pd.read_csv(path)
    if 'gene' not in df.columns:
        raise ValueError("External gene list must have a 'gene' column.")
    return df['gene'].dropna().unique().tolist()

def run_dual_go_enrichment(gene_df, out_base):
    gp = GProfiler(return_dataframe=True)
    all_genes = gene_df['name'].unique().tolist()

    if background_mode == "custom_list":
        print(f"Loading custom background gene list from: {EXTERNAL_GENE_LIST_PATH}")
        all_genes = load_external_gene_list(EXTERNAL_GENE_LIST_PATH)
        print(f"Loaded {len(all_genes)} genes from custom list.")

    for cluster_id in gene_df[clustering_type].unique():
        cluster_genes = gene_df[gene_df[clustering_type] == cluster_id]['name'].tolist()
        if not cluster_genes: continue

        for mode, bg in [('human', None), ('custom', all_genes)]:
            print(f"[Cluster {cluster_id}] GO enrichment with {mode} background...")
            domain_scope = 'custom' if mode == 'custom' else 'annotated'
            res = gp.profile(
                organism='hsapiens', query=cluster_genes, domain_scope=domain_scope,
                background=bg, user_threshold=0.05, no_evidences=False,
                significance_threshold_method='fdr'
            )
            if res.empty:
                print(f"No enrichment for cluster {cluster_id} with {mode} background.")
                continue

            res['Fold Enrichment'] = (res['intersection_size']/res['query_size']) / (res['term_size']/res['effective_domain_size'])
            res['%'] = res['intersection_size'] / res['query_size'] * 100
            res['Genes'] = res['intersections'].apply(parse_genes)
            res['Term'] = res['native'] + '~' + res['name']
            res['Category'] = res['source'].map(SOURCE_CATEGORY_MAP).fillna(res['source'])

            res.rename(columns={
                'name': 'description',
                'intersection_size': 'Count',
                'query_size': 'List Total',
                'term_size': 'Pop Total',
                'p_value': 'PValue'
            }, inplace=True)

            res['FDR'] = res['PValue']  # already adjusted

            cols = ['Term', 'description', 'Category', 'Count', '%', 'PValue', 'Genes',
                    'List Total', 'Count', 'Pop Total', 'Fold Enrichment', 'FDR']

            bg_dir = os.path.join(out_base, f"{mode}_background")
            os.makedirs(bg_dir, exist_ok=True)
            res.to_csv(os.path.join(bg_dir, f"cluster_{cluster_id}_go_enrichment.csv"), columns=cols, index=False)
            print(f"Saved enrichment results to: {bg_dir}")

#----------------------------------------------------------------------------------------------------------------#
# PLOT CIRCULAR GRAPHS PER CLUSTER
#----------------------------------------------------------------------------------------------------------------#

def save_go_legend(df, go_terms, colors, legend_path):
    info = []
    for term in go_terms:
        row = df[df['GO Term Name'] == term].iloc[0]
        cat = row['Category'].replace('GOTERM_', '').replace('_DIRECT', '')
        desc = row['Term'].split('~')[1] if '~' in row['Term'] else ''
        desc = desc[:37] + '...' if len(desc) > 40 else desc
        info.append((term, cat, desc))

    n = len(go_terms)
    row_h = 0.35
    fig_h = max(1, row_h * n + 0.5)
    fig, ax = plt.subplots(figsize=(6, fig_h))
    ax.axis('off')

    for i, (term, cat, desc) in enumerate(info):
        y = fig_h - (i + 1) * row_h - 0.1
        ax.add_patch(mpatches.Rectangle((0.05, y), 0.25, 0.2, facecolor=colors[i], edgecolor='black'))
        ax.text(0.33, y + 0.15, term, fontsize=7, va='center', ha='left', family='monospace')
        ax.text(1.1, y + 0.15, cat, fontsize=7, va='center', ha='left')
        ax.text(1.6, y + 0.15, desc, fontsize=6, va='center', ha='left', wrap=True)

    ax.set_xlim(0, 5)
    ax.set_ylim(0, fig_h)
    plt.tight_layout(pad=0.1)
    plt.savefig(legend_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved GO term legend to: {legend_path}")

def plot_go_clusters(bg_dirs):
    for bg_dir in bg_dirs:
        for file in os.listdir(bg_dir):
            if not file.endswith(".csv"):
                continue

            df = pd.read_csv(os.path.join(bg_dir, file))
            df.columns = ['Term', 'description1', 'description2', 'Category', 'Count', 'Percent',
                          'PValue', 'Genes', 'List Total', 'Count2', 'Pop Total', 'Fold Enrichment', 'FDR']

            df = df[df['Category'].isin(GO_CATEGORIES)].copy()
            df['GO ID'] = df['Term'].apply(lambda x: x.split('~')[0] if '~' in x else x)
            df['GO Term Name'] = df['GO ID']
            df['Fold Enrichment'] = pd.to_numeric(df['Fold Enrichment'], errors='coerce')
            df['FDR'] = pd.to_numeric(df['FDR'], errors='coerce')

            if len(df) > 10:
                df = df.nsmallest(10, 'FDR')

            gene_go_connections = {}
            for _, row in df.iterrows():
                genes = [g.strip() for g in row['Genes'].split(',')]
                go_term = row['GO Term Name']
                for gene in genes:
                    gene_go_connections.setdefault(gene, []).append({
                        'go_term': go_term,
                        'fold_enrichment': row['Fold Enrichment'],
                        'fdr': row['FDR']
                    })

            gene_fe = {gene: min(conns, key=lambda x: x['fdr'])['fold_enrichment']
                       for gene, conns in gene_go_connections.items()}
            gene_pvalues = {gene: min(conns, key=lambda x: x['fdr'])['fdr']
                            for gene, conns in gene_go_connections.items()}

            all_genes = sorted(gene_fe.keys(), key=lambda g: gene_fe[g])
            all_go_terms = sorted(df['GO Term Name'].unique())
            go_term_counts = df.groupby('GO Term Name')['Count'].sum().to_dict()

            radius, arc_width = 1.0, 0.1
            vmin, vmax = (min(gene_fe.values()) if gene_fe else 0, max(gene_fe.values()) if gene_fe else 1)
            cmap, norm = plt.cm.coolwarm, plt.Normalize(vmin=vmin, vmax=vmax)
            palette = [c for c in sns.color_palette("pastel", len(all_go_terms) + 3) if np.mean(c) < 0.9][:len(all_go_terms)]
            go_term_color_map = {term: palette[i] for i, term in enumerate(all_go_terms)}

            fig, ax = plt.subplots(figsize=(12, 12), subplot_kw=dict(aspect="equal"))
            ax.set(xlim=[-1.4, 1.4], ylim=[-1.4, 1.4])
            ax.axis('off')

            total_angle = np.pi
            go_scale = total_angle / sum(go_term_counts.values()) if go_term_counts else 1
            gene_scale = total_angle / len(all_genes) if all_genes else 1

            node_angles, arc_props = {}, {}
            angle = 0
            for term in all_go_terms:
                arc = go_term_counts.get(term, 0) * go_scale
                node_angles[term] = angle + arc / 2
                arc_props[term] = (angle, angle + arc, go_term_color_map[term])
                angle += arc

            angle = np.pi
            for gene in all_genes:
                arc = gene_scale
                node_angles[gene] = angle + arc / 2
                arc_props[gene] = (angle, angle + arc, cmap(norm(gene_fe[gene])))
                angle += arc

            for node, (start, end, color) in arc_props.items():
                ax.add_patch(mpatches.Wedge((0, 0), radius, np.degrees(start), np.degrees(end),
                                            width=arc_width, facecolor=color, edgecolor='black', lw=0.5))

                mid = node_angles[node]
                x, y = (radius + 0.1) * np.cos(mid), (radius + 0.1) * np.sin(mid)
                deg = np.degrees(mid)
                ha = 'right' if 90 <= deg <= 270 else 'left'
                rot = deg + 180 if ha == 'right' else deg

                label_color = 'black'
                if node in all_genes:
                    fdr = gene_pvalues.get(node, 1)
                    label_color = 'darkgreen' if fdr < 0.01 else ('orange' if fdr < 0.05 else 'purple')

                ax.text(x, y, node, ha=ha, va='center', rotation=rot, rotation_mode='anchor', fontsize=8, color=label_color)

                if node in go_term_counts:
                    x2, y2 = (radius - arc_width / 2) * np.cos(mid), (radius - arc_width / 2) * np.sin(mid)
                    ax.text(x2, y2, str(go_term_counts[node]), ha='center', va='center', fontsize=7, color='black', weight='bold')

            for _, row in df.iterrows():
                term = row['GO Term Name']
                genes = [g.strip() for g in row['Genes'].split(',')]
                x0, y0 = radius * np.cos(node_angles[term]), radius * np.sin(node_angles[term])
                color = go_term_color_map[term]
                for gene in genes:
                    if gene in node_angles:
                        x1, y1 = radius * np.cos(node_angles[gene]), radius * np.sin(node_angles[gene])
                        verts = [[x0, y0], [0, 0], [x1, y1]]
                        path = mpatches.Path(verts, [mpatches.Path.MOVETO, mpatches.Path.CURVE3, mpatches.Path.CURVE3])
                        ax.add_patch(mpatches.PathPatch(path, edgecolor=color, facecolor='none', lw=2, alpha=0.5))

            cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical', pad=0.15, shrink=0.6)
            cbar.set_label('Fold Enrichment')

            cluster_id = os.path.splitext(file)[0]
            plot_dir = os.path.join(bg_dir, 'plots')
            os.makedirs(plot_dir, exist_ok=True)
            plot_path = os.path.join(plot_dir, f"{cluster_id}_circular_plot.png")

            try:
                fig.tight_layout()
            except Exception as e:
                print(f"⚠️ Warning: tight_layout failed for cluster {cluster_id}: {e}")

            plt.savefig(plot_path, dpi=300)
            plt.close()
            print(f"✅ Saved plot to: {plot_path}")

            # Save legend if GO terms exist
            legend_path = os.path.join(plot_dir, f"{cluster_id}_go_term_legend.png")
            if all_go_terms:
                save_go_legend(df, all_go_terms, [go_term_color_map[t] for t in all_go_terms], legend_path)


# # === Main ===
# if __name__ == "__main__":
#     gene_df = load_gene_df(CSV_PATH)
#     run_dual_go_enrichment(gene_df, OUTPUT_BASE)
#     plot_go_clusters(BACKGROUND_DIRS)

#----------------------------------------------------------------------------------------------------------------#
# PROTEIN GROUPS ENRICHMENT ANALYSIS
#----------------------------------------------------------------------------------------------------------------#

# --- Constants ---
MAIN_CAT_PATH = os.path.join(DATASET_DIR, "Protein_Groups.xlsx") 
UNIPROT_MAP_PATH = os.path.join(DATASET_DIR, "uniprotkb_HUMAN_REVIEWED_GENES_REFSEQ_SEQUENCES_2025_07_01.tsv") 

# --- Load and build gene_df with Main Category ---
def load_gene_df_2(csv_path=CSV_PATH):
    cluster_file = pd.read_csv(csv_path)
    main_cat_map = pd.read_excel(MAIN_CAT_PATH)
    uniprot_map = pd.read_csv(UNIPROT_MAP_PATH, sep="\t")

    # Build gene -> UID map (first match only)
    gene_to_uid = {}
    for _, row in uniprot_map.iterrows():
        uid = row["Entry"]
        gene_names = str(row["Gene Names"]).split()
        for gene in gene_names:
            if gene not in gene_to_uid:
                gene_to_uid[gene] = uid

    # Build gene -> Main Category map
    gene_to_cat = dict(zip(main_cat_map["Gene name"], main_cat_map["Main Category"]))

    # Rename and map columns
    cluster_file.rename(columns={"name": "gene"}, inplace=True)
    cluster_file["uid"] = cluster_file["gene"].map(gene_to_uid)
    cluster_file["Main Category"] = cluster_file["gene"].map(gene_to_cat)
    cluster_file["details"] = None
    cluster_file["API_check"] = 0

    # Inference logic
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

    # Infer missing Main Categories using UniProt API
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

    return cluster_file

# --- Category enrichment function ---
def category_enrichment_analysis(gene_df, background_dir):
    #for bg_dir in background_dirs:
    output_dir = os.path.join(background_dir, "protein_group__enrichment_results")
    os.makedirs(output_dir, exist_ok=True)

    background = gene_df.dropna(subset=["Main Category"])
    M = background["gene"].nunique()
    background_cat_counts = background["Main Category"].value_counts().to_dict()

    for cluster_id, subdf in gene_df.dropna(subset=[clustering_type]).groupby(clustering_type):
        subdf = subdf.dropna(subset=["Main Category"])
        N = len(subdf["gene"].unique())
        cluster_cat_counts = subdf["Main Category"].value_counts()

        enrichment_results = []
        for cat, x in cluster_cat_counts.items():
            K = background_cat_counts.get(cat, 0)
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

        df_enrich = pd.DataFrame(enrichment_results).sort_values("PValue").reset_index(drop=True)
        m = len(df_enrich)
        df_enrich["FDR"] = df_enrich["PValue"] * m / (df_enrich.index + 1)
        df_enrich["FDR"] = df_enrich["FDR"].clip(upper=1.0)

        safe_cluster_id = str(cluster_id).replace(" ", "_")
        df_enrich.to_csv(os.path.join(output_dir, f"cluster_{safe_cluster_id}_enrichment.csv"), index=False)

        #circular_plot_path = os.path.join(output_dir, f"cluster_{safe_cluster_id}_circular_plot.png")
        create_circular_plot(output_dir)

    print(f"✅ Category enrichment complete. Results saved in {output_dir}")


def create_circular_plot(bg_dir):
    #for bg_dir in bg_dirs:
    for file in os.listdir(bg_dir):
        if not file.endswith(".csv"):
            continue

        df = pd.read_csv(os.path.join(bg_dir, file))
        df.columns = ['Term', 'description', 'Category', 'Count', '%',
                        'PValue', 'Genes', 'List Total', 'Pop Total', 'Fold Enrichment', 'FDR']

        #df = df[df['Category'].isin(GO_CATEGORIES)].copy()
        df['GO ID'] = df['Term'].apply(lambda x: x.split('~')[0] if '~' in x else x)
        df['GO Term Name'] = df['GO ID']
        df['Fold Enrichment'] = pd.to_numeric(df['Fold Enrichment'], errors='coerce')
        df['FDR'] = pd.to_numeric(df['FDR'], errors='coerce')

        if len(df) > 10:
            df = df.nsmallest(10, 'FDR')

        gene_go_connections = {}
        for _, row in df.iterrows():
            genes = [g.strip() for g in row['Genes'].split(';')]
            go_term = row['GO Term Name']
            for gene in genes:
                gene_go_connections.setdefault(gene, []).append({
                    'go_term': go_term,
                    'fold_enrichment': row['Fold Enrichment'],
                    'fdr': row['FDR']
                })

        gene_fe = {gene: min(conns, key=lambda x: x['fdr'])['fold_enrichment']
                    for gene, conns in gene_go_connections.items()}
        gene_fdr = {gene: min(conns, key=lambda x: x['fdr'])['fdr']
                    for gene, conns in gene_go_connections.items()}

        all_genes = sorted(gene_fe.keys(), key=lambda g: gene_fe[g])
        all_go_terms = sorted(df['GO Term Name'].unique())
        go_term_counts = df.groupby('GO Term Name')['Count'].sum().to_dict()

        radius = 1.0
        arc_width = 0.1

        vmin = min(gene_fe.values()) if gene_fe else 0
        vmax = max(gene_fe.values()) if gene_fe else 1
        cmap = plt.cm.coolwarm
        norm = plt.Normalize(vmin=vmin, vmax=vmax)

        palette = [c for c in sns.color_palette("pastel", len(all_go_terms) + 3) if np.mean(c) < 0.9][:len(all_go_terms)]
        go_term_color_map = {term: palette[i] for i, term in enumerate(all_go_terms)}

        fig, ax = plt.subplots(figsize=(12, 12), subplot_kw=dict(aspect="equal"))
        ax.set_xlim([-1.4, 1.4])
        ax.set_ylim([-1.4, 1.4])
        ax.axis('off')

        total_angle = np.pi
        go_scale = total_angle / sum(go_term_counts.values()) if go_term_counts else 1
        gene_scale = total_angle / len(all_genes) if all_genes else 1

        node_angles = {}
        arc_props = {}

        angle = 0
        for term in all_go_terms:
            arc = go_term_counts.get(term, 0) * go_scale
            node_angles[term] = angle + arc / 2
            arc_props[term] = (angle, angle + arc, go_term_color_map[term])
            angle += arc

        angle = np.pi
        for gene in all_genes:
            arc = gene_scale
            node_angles[gene] = angle + arc / 2
            arc_props[gene] = (angle, angle + arc, cmap(norm(gene_fe[gene])))
            angle += arc

        for node, (start, end, color) in arc_props.items():
            ax.add_patch(mpatches.Wedge((0, 0), radius, np.degrees(start), np.degrees(end),
                                        width=arc_width, facecolor=color, edgecolor='black', lw=0.5))

            mid = node_angles[node]
            x, y = (radius + 0.1) * np.cos(mid), (radius + 0.1) * np.sin(mid)
            deg = np.degrees(mid)
            ha = 'right' if 90 <= deg <= 270 else 'left'
            rot = deg + 180 if ha == 'right' else deg

            if node in all_genes:
                fdr = gene_fdr.get(node, 1)
                label_color = 'darkgreen' if fdr < 0.01 else ('orange' if fdr < 0.05 else 'purple')
            else:
                label_color = 'black'

            ax.text(x, y, node, ha=ha, va='center', rotation=rot, rotation_mode='anchor', fontsize=8, color=label_color)

            if node in go_term_counts:
                x2, y2 = (radius - arc_width / 2) * np.cos(mid), (radius - arc_width / 2) * np.sin(mid)
                ax.text(x2, y2, str(go_term_counts[node]), ha='center', va='center', fontsize=7, color='black', weight='bold')

        for _, row in df.iterrows():
            term = row['GO Term Name']
            genes = [g.strip() for g in row['Genes'].split(';')]
            x0, y0 = radius * np.cos(node_angles[term]), radius * np.sin(node_angles[term])
            color = go_term_color_map[term]

            for gene in genes:
                if gene in node_angles:
                    x1, y1 = radius * np.cos(node_angles[gene]), radius * np.sin(node_angles[gene])
                    verts = [[x0, y0], [0, 0], [x1, y1]]
                    path = mpatches.Path(verts, [mpatches.Path.MOVETO, mpatches.Path.CURVE3, mpatches.Path.CURVE3])
                    ax.add_patch(mpatches.PathPatch(path, edgecolor=color, facecolor='none', lw=2, alpha=0.5))

        cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical', pad=0.15, shrink=0.6)
        cbar.set_label('Fold Enrichment')

        cluster_id = os.path.splitext(file)[0]
        plot_dir = os.path.join(bg_dir, 'plots')
        os.makedirs(plot_dir, exist_ok=True)

        plot_path = os.path.join(plot_dir, f"{cluster_id}_circular_plot.png")
        plt.tight_layout()
        plt.savefig(plot_path, dpi=300)
        plt.close()
        print(f"Saved plot to: {plot_path}")

        legend_path = os.path.join(plot_dir, f"{cluster_id}_go_term_legend.png")
        save_go_legend(df, all_go_terms, [go_term_color_map[t] for t in all_go_terms], legend_path)

#----------------------------------------------------------------------------------------------------------------#
# CORUM DB COMPLEXES ENRICHMENT ANALYSIS
#----------------------------------------------------------------------------------------------------------------#

CORUM_CSV_PATH= os.path.join(DATASET_DIR, "corum_humanComplexes_7_2025.csv") 

def harmonic_mean(a, b):
    return 2 * a * b / (a + b) if (a + b) > 0 else 0

def precision_recall_product(p, b):
    return (len(p & b) ** 2) / (len(p) * len(b)) if p and b else 0

def jaccard(p, b):
    return len(p & b) / len(p | b) if p | b else 0

def corum_enrichment_analysis(cluster_file_path, corum_csv_path, output_dir):
    df = pd.read_csv(cluster_file_path)
    df = df.dropna(subset=[clustering_type, "name"])
    df["genes"] = df["name"].astype(str).str.split(";|,")
    df = df.explode("genes")
    df["genes"] = df["genes"].str.strip()

    clusters = df.groupby(clustering_type)["genes"].apply(set).to_dict()

    corum_df = pd.read_csv(corum_csv_path, encoding="ISO-8859-1")
    corum_df = corum_df[corum_df["organism"] == "Human"]
    corum_df["genes"] = corum_df["subunits_gene_name"].fillna("").astype(str).str.split(";").apply(lambda x: set(map(str.strip, x)))
    corum_complexes = corum_df.to_dict("records")

    all_genes = set(df["genes"].dropna())
    M = len(all_genes)

    pvals = []
    enrichment_results = []

    for cluster_id, cluster_genes in clusters.items():
        cluster_results = []
        for complex_row in corum_complexes:
            complex_genes = complex_row["genes"]
            if not complex_genes:
                continue
            overlap = cluster_genes & complex_genes
            if not overlap:
                continue

            N = len(cluster_genes)
            K = len(complex_genes)
            x = len(overlap)
            pval = hypergeom.sf(x - 1, M, K, N)
            pvals.append(pval)

            precision = x / N
            recall = x / K
            fscore = harmonic_mean(precision, recall)
            jac = jaccard(cluster_genes, complex_genes)

            cluster_results.append({
                "Term": complex_row["complex_name"],
                "description": complex_row["complex_name"],
                "Category": "complex_name",
                "Count": x,
                "%": round(x / N * 100, 2),
                "PValue": pval,
                "Genes": "; ".join(overlap),
                "List Total": N,
                "Pop Total": M,
                "Fold Enrichment": round((x / N) / (K / M), 2) if K > 0 else 0,
                "Jaccard": round(jac, 3),
                "F-measure": round(fscore, 3)
            })

        # FDR correction
        if cluster_results:
            pvals_cluster = [res["PValue"] for res in cluster_results]
            _, fdrs, _, _ = multipletests(pvals_cluster, method="fdr_bh")
            for res, fdr in zip(cluster_results, fdrs):
                res["FDR"] = round(min(fdr, 1.0), 5)
            df_out = pd.DataFrame(cluster_results).sort_values("FDR").reset_index(drop=True)
            safe_cluster_id = str(cluster_id).replace(" ", "_")
            cluster_output_path = os.path.join(output_dir, f"cluster_{safe_cluster_id}_enrichment.csv")
            df_out.to_csv(cluster_output_path, index=False)
            print(f"✅ CORUM enrichment for cluster {cluster_id} saved to {cluster_output_path}")


def save_go_legend(df, all_go_terms, colors, legend_path):
    fig, ax = plt.subplots(figsize=(3, len(all_go_terms) * 0.3))
    for i, (term, color) in enumerate(zip(all_go_terms, colors)):
        ax.barh(i, 1, color=color)
        ax.text(1.05, i, term, va='center')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(legend_path)
    plt.close()

def create_circular_plot_top10_CORUM(enrich_dir):
    # Find cluster files and their minimal FDR
    cluster_fdrs = []
    cluster_files = []

    for file in os.listdir(enrich_dir):
        if not file.endswith(".csv"):
            continue
        df = pd.read_csv(os.path.join(enrich_dir, file))
        if 'FDR' not in df.columns:
            continue
        df['FDR'] = pd.to_numeric(df['FDR'], errors='coerce')
        min_fdr = df['FDR'].min()
        if pd.notnull(min_fdr):
            cluster_fdrs.append(min_fdr)
            cluster_files.append(file)

    # Select top 10 clusters by minimal FDR
    top10_indices = np.argsort(cluster_fdrs)[:10]
    top10_files = [cluster_files[i] for i in top10_indices]

    for file in top10_files:
        df = pd.read_csv(os.path.join(enrich_dir, file))
        df.columns = ['Term', 'description', 'Category', 'Count', '%',
                      'PValue', 'Genes', 'List Total', 'Pop Total', 'Fold Enrichment', 'Jaccard', 'F-measure', 'FDR']

        # Keep only top 10 terms by FDR per cluster (optional)
        if len(df) > 10:
            df = df.nsmallest(10, 'FDR')

        df['GO Term Name'] = df['Term']
        df['Fold Enrichment'] = pd.to_numeric(df['Fold Enrichment'], errors='coerce')
        df['FDR'] = pd.to_numeric(df['FDR'], errors='coerce')

        gene_go_connections = {}
        for _, row in df.iterrows():
            genes = [g.strip() for g in row['Genes'].split(';')]
            go_term = row['GO Term Name']
            for gene in genes:
                gene_go_connections.setdefault(gene, []).append({
                    'go_term': go_term,
                    'fold_enrichment': row['Fold Enrichment'],
                    'fdr': row['FDR']
                })

        gene_fe = {gene: min(conns, key=lambda x: x['fdr'])['fold_enrichment']
                    for gene, conns in gene_go_connections.items()}
        gene_fdr = {gene: min(conns, key=lambda x: x['fdr'])['fdr']
                    for gene, conns in gene_go_connections.items()}

        all_genes = sorted(gene_fe.keys(), key=lambda g: gene_fe[g])
        all_go_terms = sorted(df['GO Term Name'].unique())
        go_term_counts = df.groupby('GO Term Name')['Count'].sum().to_dict()

        radius = 1.0
        arc_width = 0.1

        vmin = min(gene_fe.values()) if gene_fe else 0
        vmax = max(gene_fe.values()) if gene_fe else 1
        cmap = plt.cm.coolwarm
        norm = plt.Normalize(vmin=vmin, vmax=vmax)

        palette = [c for c in sns.color_palette("pastel", len(all_go_terms) + 3) if np.mean(c) < 0.9][:len(all_go_terms)]
        go_term_color_map = {term: palette[i] for i, term in enumerate(all_go_terms)}

        fig, ax = plt.subplots(figsize=(12, 12), subplot_kw=dict(aspect="equal"))
        ax.set_xlim([-1.4, 1.4])
        ax.set_ylim([-1.4, 1.4])
        ax.axis('off')

        total_angle = np.pi
        go_scale = total_angle / sum(go_term_counts.values()) if go_term_counts else 1
        gene_scale = total_angle / len(all_genes) if all_genes else 1

        node_angles = {}
        arc_props = {}

        angle = 0
        for term in all_go_terms:
            arc = go_term_counts.get(term, 0) * go_scale
            node_angles[term] = angle + arc / 2
            arc_props[term] = (angle, angle + arc, go_term_color_map[term])
            angle += arc

        angle = np.pi
        for gene in all_genes:
            arc = gene_scale
            node_angles[gene] = angle + arc / 2
            arc_props[gene] = (angle, angle + arc, cmap(norm(gene_fe[gene])))
            angle += arc

        for node, (start, end, color) in arc_props.items():
            ax.add_patch(mpatches.Wedge((0, 0), radius, np.degrees(start), np.degrees(end),
                                        width=arc_width, facecolor=color, edgecolor='black', lw=0.5))

            mid = node_angles[node]
            x, y = (radius + 0.1) * np.cos(mid), (radius + 0.1) * np.sin(mid)
            deg = np.degrees(mid)
            ha = 'right' if 90 <= deg <= 270 else 'left'
            rot = deg + 180 if ha == 'right' else deg

            if node in all_genes:
                fdr = gene_fdr.get(node, 1)
                label_color = 'darkgreen' if fdr < 0.01 else ('orange' if fdr < 0.05 else 'purple')
            else:
                label_color = 'black'

            ax.text(x, y, node, ha=ha, va='center', rotation=rot, rotation_mode='anchor', fontsize=8, color=label_color)

            if node in go_term_counts:
                x2, y2 = (radius - arc_width / 2) * np.cos(mid), (radius - arc_width / 2) * np.sin(mid)
                ax.text(x2, y2, str(go_term_counts[node]), ha='center', va='center', fontsize=7, color='black', weight='bold')

        for _, row in df.iterrows():
            term = row['GO Term Name']
            genes = [g.strip() for g in row['Genes'].split(';')]
            x0, y0 = radius * np.cos(node_angles[term]), radius * np.sin(node_angles[term])
            color = go_term_color_map[term]

            for gene in genes:
                if gene in node_angles:
                    x1, y1 = radius * np.cos(node_angles[gene]), radius * np.sin(node_angles[gene])
                    verts = [[x0, y0], [0, 0], [x1, y1]]
                    path = mpatches.Path(verts, [mpatches.Path.MOVETO, mpatches.Path.CURVE3, mpatches.Path.CURVE3])
                    ax.add_patch(mpatches.PathPatch(path, edgecolor=color, facecolor='none', lw=2, alpha=0.5))

        cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical', pad=0.15, shrink=0.6)
        cbar.set_label('Fold Enrichment')

        cluster_id = os.path.splitext(file)[0]
        plot_dir = os.path.join(enrich_dir, 'plots')
        os.makedirs(plot_dir, exist_ok=True)

        plot_path = os.path.join(plot_dir, f"{cluster_id}_circular_plot.png")
        plt.tight_layout()
        plt.savefig(plot_path, dpi=300)
        plt.close()
        print(f"Saved plot to: {plot_path}")

        legend_path = os.path.join(plot_dir, f"{cluster_id}_go_term_legend.png")
        save_go_legend(df, all_go_terms, [go_term_color_map[t] for t in all_go_terms], legend_path)

if __name__ == "__main__":

    if group_enrichment_mode == 0:
        gene_df = load_gene_df(CSV_PATH)
        run_dual_go_enrichment(gene_df, OUTPUT_BASE)
        plot_go_clusters(BACKGROUND_DIRS)
    elif group_enrichment_mode == 1:
        gene_df = load_gene_df(CSV_PATH)
        run_dual_go_enrichment(gene_df, OUTPUT_BASE)
        plot_go_clusters(BACKGROUND_DIRS)
        gene_df_2 = load_gene_df_2(CSV_PATH)
        category_enrichment_analysis(gene_df_2, OUTPUT_BASE)
        output_dir_corum = os.path.join(OUTPUT_BASE, "corum_enrichment_results")
        os.makedirs(output_dir_corum, exist_ok=True)
        corum_enrichment_analysis(CSV_PATH, CORUM_CSV_PATH, output_dir_corum)
        create_circular_plot_top10_CORUM(output_dir_corum)
    
    # elif group_enrichment_mode == 2:
    #     output_dir = os.path.join(OUTPUT_BASE, "corum_enrichment_results")
    #     os.makedirs(output_dir, exist_ok=True)
    #     corum_enrichment_analysis(CSV_PATH, CORUM_CSV_PATH, output_dir)
    #     create_circular_plot_top10_CORUM(output_dir)

        
    else:
        raise ValueError("Invalid mode selected. Use 0 for original workflow or 1 for category enrichment.")
