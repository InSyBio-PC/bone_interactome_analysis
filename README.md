# Bone Interactome Analysis

Computational reconstruction and analysis of a bone degeneration-related protein interactome using machine-learning-based protein–protein interaction prediction, network reconstruction, functional analysis, and molecular dynamics.

> [!IMPORTANT]
> This repository is organized into two major parts:
>
> **Part 1 — Bone EOA Analysis:** reconstruction and systems-level analysis of the bone protein–protein interaction network.
>
> **Part 2 — Molecular Dynamics Analysis:** structural and molecular-dynamics investigation of selected interactions from the reconstructed network.
>
> Large datasets and Cytoscape project files that cannot be hosted directly on GitHub are available separately. See [Large Files and Datasets](#large-files-and-datasets).

---

# Contents

- [Part 1 — Bone EOA Analysis](#part-1--bone-eoa-analysis)
  - [Overview](#part-1-overview)
  - [1. Initial Bone Protein Dataset](#11-initial-bone-protein-dataset)
  - [2. Candidate Bone PPI Generation](#12-candidate-bone-ppi-generation)
  - [3. PPI Feature Calculation](#13-ppi-feature-calculation)
  - [4. EOA Classification and Regression](#14-eoa-classification-and-regression)
  - [5. Positive PPI Filtering](#15-positive-ppi-filtering)
  - [6. Interactor Expansion](#16-interactor-expansion)
  - [7. Interaction Mapping](#17-interaction-mapping)
  - [8. Final Bone PPI Network](#18-final-bone-ppi-network)
  - [9. Network Reconstruction and Clustering](#19-network-reconstruction-and-clustering)
  - [10. Network Analytics](#110-network-analytics)
  - [11. Functional Enrichment Analysis](#111-functional-enrichment-analysis)
  - [Part 1 Dataset Summary](#part-1-dataset-summary)
  - [Part 1 Repository Structure](#part-1-repository-structure)
- [Part 2 — Molecular Dynamics Analysis](#part-2--molecular-dynamics-analysis)
- [Large Files and Datasets](#large-files-and-datasets)
- [Reproducing the Analysis](#reproducing-the-analysis)
- [Acknowledgements](#acknowledgements)

---

# Part 1 — Bone EOA Analysis

<a id="part-1-overview"></a>

## Overview

Part 1 reconstructs a bone degeneration-related protein–protein interaction network starting from a curated set of **44 osteoporosis-related proteins (OPs)**.

The workflow combines:

- generation of candidate binary PPIs;
- calculation of PPI features;
- EOA-based classification and regression;
- confidence-based interaction filtering;
- expansion to the interactors of the original OP network;
- STRING DB and iRefIndex mapping;
- gene and OP annotation;
- final PPI network reconstruction;
- Cytoscape-based clustering;
- network topology analysis; and
- functional enrichment analysis.

### Part 1 at a glance

| Stage | Result |
|---|---:|
| Osteoporosis-related proteins | **44** |
| Reviewed human UniProt entries | **20,421** |
| Candidate OP PPIs | **877,501** |
| EOA-positive OP PPIs | **479,022** |
| Filtered OP PPIs | **1,079** |
| OP interactors | **995** |
| Candidate OP-interactor combinations | **494,515** |
| EOA-positive interactor PPIs | **46,722** |
| Filtered interactor PPIs | **3,284** |
| **Final network** | **4,363 PPIs** |

---

## 1.1 Initial Bone Protein Dataset

[KEEP THE TEXT FROM THE PREVIOUS README SECTION HERE]

---

## 1.2 Candidate Bone PPI Generation

[KEEP THE TEXT FROM THE PREVIOUS README SECTION HERE]

---

## 1.3 PPI Feature Calculation

[KEEP THE TEXT FROM THE PREVIOUS README SECTION HERE]

The feature-calculation methodology used in this analysis follows the workflow implemented in:

**[TR_PPI — PPI Feature Calculation Workflow](https://github.com/HarrisZavs/TR_PPI/tree/main)**

---

## 1.4 EOA Classification and Regression

[KEEP THE CLASSIFICATION/REGRESSION SECTION FROM THE PREVIOUS VERSION]

---

## 1.5 Positive PPI Filtering

[KEEP THE FILTERING SECTION AND TABLES FROM THE PREVIOUS VERSION]

> [!IMPORTANT]
> **Both filtering criteria were required:**
>
> `Probability Score > 0.517074`
>
> **AND**
>
> `mean_prob_aff > 0.459257`

---

## 1.6 Interactor Expansion

[KEEP THE INTERACTOR EXPANSION SECTION]

---

## 1.7 Interaction Mapping

[KEEP THE STRING DB / iRefIndex / UniProt MAPPING SECTION]

---

## 1.8 Final Bone PPI Network

[KEEP THE FINAL 4,363 PPI NETWORK SECTION]

---

## 1.9 Network Reconstruction and Clustering

[KEEP THE CYTOSCAPE + MCL / MCODE / LEIDEN / GLAY SECTION AND TABLE]

---

## 1.10 Network Analytics

Network topology and clustering properties were assessed using Cytoscape-exported network attributes and Python-based analysis.

The calculated properties include:

- density;
- degree;
- clustering coefficient;
- betweenness centrality;
- closeness centrality;
- cluster density;
- cluster inertia; and
- silhouette-based cluster separation.

The corresponding files are located under:

```text
bone_EOA_analysis/network_reconstruction/network_analytics/
```

---

## 1.11 Functional Enrichment Analysis

Functional enrichment analysis was performed on clusters derived from the selected **MCL granularity-4** network.

The analyses include:

- Gene Ontology enrichment;
- Biological Process enrichment;
- Molecular Function enrichment;
- Cellular Component enrichment;
- protein-group enrichment; and
- CORUM-related enrichment.

Where applicable, statistically significant enrichment was defined as:

```text
Adjusted P-value < 0.05
```

Results, plots, and analysis scripts are located under:

```text
bone_EOA_analysis/network_reconstruction/MCL_gran_4_cluster_enrichment/
```

---

## Part 1 Dataset Summary

| Dataset | Rows | Description |
|---|---:|---|
| `bone_PPI_combs_raw.csv` | 877,501 | Raw features for candidate OP PPIs |
| `bone_PPI_combs_processed.csv` | 877,501 | Processed features for candidate OP PPIs |
| `bone_interactor_combs_raw.csv` | 494,515 | Raw features for candidate interactor PPIs |
| `bone_interactor_combs_processed.csv` | 494,515 | Processed features for candidate interactor PPIs |
| `bone_ppis_EOA_all_predictions_mapped.csv` | 877,501 | Complete mapped OP PPI predictions |
| `bone_interactor_combs_EOA_all_predictions_mapped.csv` | 494,515 | Complete mapped interactor predictions |
| `BOTH_FILTER_bone_ppis_EOA_positive_fully_mapped.csv` | 1,079 | Final filtered OP PPIs |
| `BOTH_FILTER_bone_interactor_combs_EOA_positive_fully_mapped.csv` | 3,284 | Final filtered interactor PPIs |
| `final_bone_ppi_network_file.csv` | **4,363** | **Final network used for reconstruction** |

---

## Part 1 Repository Structure

```text
bone_EOA_analysis/
│
├── Datasets/
│   ├── reference/
│   ├── feature_calculated/
│   ├── intermediate/
│   └── results/
│
├── fc_EOA_prediction_codes/
│   ├── prediction scripts
│   └── models/
│
├── mapping_codes/
│   ├── bone_ppi_mappings.ipynb
│   └── bone_interactor_combs_mappings.ipynb
│
└── network_reconstruction/
    ├── network figures
    ├── network_analytics/
    │   ├── codes/
    │   └── network-analysis outputs
    │
    └── MCL_gran_4_cluster_enrichment/
        ├── codes/
        ├── GO enrichment results
        ├── protein-group enrichment results
        ├── CORUM enrichment results
        └── plots/
```

---

# Part 2 — Molecular Dynamics Analysis

> [!NOTE]
> **Part 2** contains the molecular dynamics component of the Bone Interactome Analysis.

---

# Large Files and Datasets

The large datasets required for the analyses, together with the Cytoscape project file, are available separately:

### [Google Drive — Bone Interactome Large Files](PASTE_GOOGLE_DRIVE_FOLDER_LINK_HERE)

> [!WARNING]
> These files are intentionally excluded from Git because of their size. They should be downloaded and restored to the corresponding repository directories when required.

The Google Drive archive should follow the repository organization:

| Google Drive folder | Contents | Corresponding repository location |
|---|---|---|
| `Part_1_Bone_EOA/feature_calculated/` | Raw and processed EOA feature matrices | `bone_EOA_analysis/Datasets/feature_calculated/` |
| `Part_1_Bone_EOA/intermediate/` | Intermediate prediction/mapping datasets | `bone_EOA_analysis/Datasets/intermediate/` |
| `Part_1_Bone_EOA/reference/` | Large reference datasets | `bone_EOA_analysis/Datasets/reference/` |
| `Part_1_Bone_EOA/results/` | Large result datasets excluded from Git | `bone_EOA_analysis/Datasets/results/` |
| `Part_1_Bone_EOA/cytoscape/` | Cytoscape project | `bone_EOA_analysis/network_reconstruction/` |
| `Part_2_Molecular_Dynamics/` | Large files associated with Part 2 | Corresponding Part 2 directories |

> [!TIP]
> Preserve the directory organization when downloading the files so that the large datasets can be restored to their expected analysis locations.

---

# Reproducing the Analysis

## Part 1

```text
44 osteoporosis-related proteins
              │
              ▼
20,421 reviewed human UniProt entries
              │
              ▼
877,501 candidate OP PPIs
              │
              ▼
PPI feature calculation
              │
              ▼
EOA classification + regression
              │
              ▼
479,022 predicted positive PPIs
              │
              ▼
Dual confidence filtering
              │
              ▼
1,079 filtered OP PPIs
              │
              ▼
995 OP interactors
              │
              ▼
494,515 candidate interactor combinations
              │
              ▼
EOA classification + regression
              │
              ▼
46,722 predicted positive PPIs
              │
              ▼
Dual confidence filtering
              │
              ▼
3,284 filtered interactor PPIs
              │
              ▼
4,363 final PPIs
              │
              ▼
Interaction and gene mapping
              │
              ▼
Network reconstruction
              │
              ▼
Clustering comparison
              │
              ▼
MCL granularity 4
              │
              ▼
Network analytics
              │
              ▼
Functional enrichment
```

## Part 2

Molecular Dynamics Analysis.

---

# Acknowledgements

The present work has been developed as part of the **REGENERATION project**, funded by the European Union’s Horizon 2020 research and innovation program under the **Marie Sklodowska-Curie RISE (Grant Agreement No. 101131255)**.

This work was supported by the **Swiss State Secretariat for Education, Research and Innovation (SERI)** under contract number **23.0086**.
