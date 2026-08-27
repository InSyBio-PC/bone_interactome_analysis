# Bone Interactome Analysis

Computational reconstruction and analysis of a **bone degeneration-related protein interactome**, combining machine-learning-based protein–protein interaction (PPI) prediction, interaction database mapping, network reconstruction, clustering, and functional enrichment analysis.

> [!IMPORTANT]
> This repository contains the code, final compact datasets, network-analysis outputs, and supporting files required to document the analysis.
>
> Several large datasets and the Cytoscape project file exceed practical GitHub repository size limits and are therefore stored separately on Google Drive. See [Large Files and Datasets](#large-files-and-datasets).

---

## Contents

- [Project Overview](#project-overview)
- [Analysis Structure](#analysis-structure)
  - [Part 1 — Bone EOA Analysis](#part-1--bone-eoa-analysis)
  - [Part 2](#part-2)
- [Part 1 Workflow](#part-1-workflow)
  - [1. Initial Bone Protein Dataset](#1-initial-bone-protein-dataset)
  - [2. Candidate Bone PPI Generation](#2-candidate-bone-ppi-generation)
  - [3. PPI Feature Calculation](#3-ppi-feature-calculation)
  - [4. EOA Classification and Regression](#4-eoa-classification-and-regression)
  - [5. Positive PPI Filtering](#5-positive-ppi-filtering)
  - [6. Interactor Expansion](#6-interactor-expansion)
  - [7. Interaction Mapping](#7-interaction-mapping)
  - [8. Final Bone PPI Network](#8-final-bone-ppi-network)
  - [9. Network Reconstruction and Clustering](#9-network-reconstruction-and-clustering)
  - [10. Network Analytics and Enrichment](#10-network-analytics-and-enrichment)
- [Key Dataset Sizes](#key-dataset-sizes)
- [Dataset Columns](#dataset-columns)
- [Repository Structure](#repository-structure)
- [Large Files and Datasets](#large-files-and-datasets)
- [Reproducing the Analysis](#reproducing-the-analysis)
- [Acknowledgements](#acknowledgements)

---

# Project Overview

The analysis starts from a curated set of **44 osteoporosis-related proteins (OPs)** and expands them against **20,421 reviewed human UniProt entries** (December 2024).

This generated:

**877,501 candidate bone PPIs**

These candidate interactions were evaluated using the **Evolutionary Optimization Algorithm (EOA)-based PPI prediction methodology**, combining classification and regression predictions.

The workflow subsequently expands the network through the predicted interactors of the OP-associated proteins and produces a final network containing:

| Network component | PPIs |
|---|---:|
| Filtered OP PPIs | 1,079 |
| Filtered OP-interactor PPIs | 3,284 |
| **Final network** | **4,363** |

The final interaction network was reconstructed and analyzed using Cytoscape and Python-based network-analysis workflows.

---

# Analysis Structure

## Part 1 — Bone EOA Analysis

Part 1 contains the complete workflow for:

- candidate PPI generation;
- PPI feature calculation;
- EOA classification and regression;
- confidence-based PPI filtering;
- OP-interactor network expansion;
- STRING DB and iRefIndex mapping;
- gene and OP annotation;
- final network generation;
- Cytoscape network reconstruction;
- network clustering;
- network analytics; and
- functional enrichment analysis.

The corresponding files are located primarily under:

```text
bone_EOA_analysis/
```

---

## Part 2

---

# Part 1 Workflow

## 1. Initial Bone Protein Dataset

A curated list of **44 OPs (osteoporosis-related proteins)** was used as the starting point of the analysis.

The corresponding reference dataset is:

```text
Datasets/reference/bone_proteins_04.xlsx
```

It contains the protein names and corresponding UniProt identifiers used to define the initial bone-related protein set.

The 44 OPs were combined with **20,421 reviewed human UniProt entries**, producing the candidate interaction space used for PPI prediction.

---

## 2. Candidate Bone PPI Generation

Combining the 44 OPs with the reviewed human UniProt proteome generated:

> **877,501 candidate bone PPIs**

The corresponding candidate-combination dataset is:

```text
Datasets/reference/bone_comb.csv
```

These candidate PPIs represent possible binary interactions involving at least one OP.

---

## 3. PPI Feature Calculation

Each candidate PPI was represented using physicochemical, sequence-derived, functional, localization, evolutionary, expression, and interaction-related features required by the EOA prediction methodology.

The feature-calculation workflow used for this project follows the methodology implemented in the **TR_PPI repository**:

[TR_PPI — PPI feature calculation workflow](https://github.com/HarrisZavs/TR_PPI/tree/main)

Feature-calculated datasets include raw and processed versions.

```text
Datasets/feature_calculated/
├── bone_PPI_combs_raw.csv
├── bone_PPI_combs_processed.csv
├── bone_interactor_combs_raw.csv
└── bone_interactor_combs_processed.csv
```

> [!NOTE]
> The `_raw` datasets contain the unprocessed calculated PPI features.
>
> The `_processed` datasets contain the processed feature representations used for prediction, including KNN imputation and arithmetic sample-wise normalization.

These files are too large for the main GitHub repository and are provided separately. See [Large Files and Datasets](#large-files-and-datasets).

---

## 4. EOA Classification and Regression

EOA prediction consists of two complementary components:

| Component | Output |
|---|---|
| Classification | Predicted interaction class and probability |
| Regression | Predicted interaction affinity |
| Combined score | Mean classifier probability / predicted affinity score |

The relevant prediction implementation is located under:

```text
bone_EOA_analysis/fc_EOA_prediction_codes/
```

The classifier determines whether a candidate pair is predicted as a positive or negative PPI.

The regression component provides a predicted interaction-affinity value.

The combined metric used downstream is:

```text
mean_prob_aff
```

representing the combined classifier probability and predicted affinity information.

---

## 5. Positive PPI Filtering

### Initial OP PPI predictions

From the original:

**877,501 candidate OP PPIs**

EOA predicted:

**479,022 positive interactions**

The classifier probability distribution among these positive predictions was:

| Statistic | Probability Score |
|---|---:|
| Count | 479,022 |
| Mean | 0.516263 |
| SD | 0.004065 |
| Minimum | 0.502439 |
| 25% | 0.517073 |
| Median | 0.517073 |
| 75% | 0.517073 |
| Maximum | 0.970732 |

The `mean_prob_aff` distribution was:

| Statistic | mean_prob_aff |
|---|---:|
| Count | 479,022 |
| Mean | 0.458876 |
| SD | 0.002052 |
| Minimum | 0.451872 |
| 25% | 0.459257 |
| Median | 0.459257 |
| 75% | 0.459290 |
| Maximum | 0.692132 |

Based on these distributions, a combined filtering strategy was applied.

> [!IMPORTANT]
> **Final PPI filtering criteria**
>
> - `Probability Score > 0.517074`
> - `mean_prob_aff > 0.459257`
>
> Both conditions must be satisfied.

Application of both filters reduced the 479,022 predicted positive PPIs to:

> **1,079 high-confidence OP PPIs**

---

## 6. Interactor Expansion

The 1,079 filtered OP PPIs contained **995 OP interactors**.

All possible combinations among the relevant OP-interactor space generated:

> **494,515 candidate OP-interactor combinations**

EOA prediction identified:

> **46,722 positive interactions**

The same dual filtering criteria were then applied:

```text
Probability Score > 0.517074
AND
mean_prob_aff > 0.459257
```

resulting in:

> **3,284 filtered positive OP-interactor PPIs**

Therefore:

```text
1,079 OP PPIs
+
3,284 OP-interactor PPIs
=
4,363 final PPIs
```

---

## 7. Interaction Mapping

EOA predictions were subsequently mapped to external interaction and annotation resources.

The mapping notebooks are located under:

```text
bone_EOA_analysis/mapping_codes/
```

Reference datasets include:

```text
Datasets/reference/
├── bone_proteins_04.xlsx
├── bone_comb.csv
├── uid_to_Gene_Names.csv
├── stringdb_ppis_curated.csv
└── irefindex_v3.csv
```

### STRING DB

`stringdb_ppis_curated.csv` contains curated human binary PPIs from **STRING DB v12.0** as used in this analysis.

Relevant fields include:

| Column | Description |
|---|---|
| `string_A` | STRING identifier of protein A |
| `string_B` | STRING identifier of protein B |
| `score` | STRING interaction score |
| `uidA` | UniProt ID of protein A |
| `uidB` | UniProt ID of protein B |
| `score_norm` | Interaction score normalized to 0–1 |

### iRefIndex

`irefindex_v3.csv` contains filtered human binary PPIs represented using UniProt identifiers and incorporates the available experimental PPI-detection methods retained during dataset preparation.

### UniProt → Gene mapping

```text
uid_to_Gene_Names.csv
```

maps UniProt identifiers to their corresponding gene names.

---

## 8. Final Bone PPI Network

The two filtered PPI populations were combined into:

```text
Datasets/results/final_bone_ppi_network_file.csv
```

containing:

> **4,363 PPIs**

The `OP_check` field distinguishes the origin of each interaction:

| `OP_check` | Meaning |
|---:|---|
| `1` | OP PPI — at least one protein belongs to the original OP set |
| `0` | OP-interactor PPI |

The two complete mapped EOA prediction datasets are:

```text
bone_ppis_EOA_all_predictions_mapped.csv
bone_interactor_combs_EOA_all_predictions_mapped.csv
```

and the filtered datasets are:

```text
BOTH_FILTER_bone_ppis_EOA_positive_fully_mapped.csv
BOTH_FILTER_bone_interactor_combs_EOA_positive_fully_mapped.csv
```

---

## 9. Network Reconstruction and Clustering

Network reconstruction and visualization were performed using **Cytoscape v3.10.1**.

Multiple clustering approaches and parameter configurations were evaluated:

| Algorithm | Configuration |
|---|---|
| MCL | Granularity 4 |
| MCL | Granularity 2.5 |
| MCL | Granularity 2 |
| MCODE | Degree coefficient 3 |
| MCODE | Degree coefficient 2 |
| Leiden | Resolution 0.1 |
| Leiden | Resolution 0.05 |
| GLay | Default |

The evaluated network metrics included:

- clustering coefficient;
- density;
- heterogeneity;
- degree; and
- silhouette coefficient.

### Clustering comparison

| Algorithm | Configuration | Clustering Coefficient | Density | Heterogeneity | Degree | Silhouette |
|---|---|---:|---:|---:|---:|---:|
| MCL | Gran. 4 | 0.263 | **0.741** | 2.025 | 2.506 | **0.586** |
| MCL | Gran. 2.5 | 0.304 | 0.706 | 2.406 | 2.899 | 0.520 |
| MCL | Gran. 2 | 0.408 | 0.491 | 2.540 | 3.980 | 0.409 |
| MCODE | DC 3 | 0.660 | 0.660 | 0.910 | 3.070 | 0.388 |
| MCODE | DC 2 | 0.270 | 0.193 | 1.599 | 3.091 | 0.345 |
| Leiden | Res. 0.1 | 0.264 | 0.504 | 2.560 | 3.141 | **0.595** |
| Leiden | Res. 0.05 | 0.310 | 0.424 | 2.677 | 3.731 | 0.553 |
| GLay | Default | 0.462 | 0.127 | 1.782 | 5.416 | 0.451 |

Leiden at resolution 0.1 achieved the highest silhouette coefficient (`0.595`), closely followed by MCL granularity 4 (`0.586`).

However, **MCL granularity 4** also achieved the highest density (`0.741`) and was therefore selected for downstream network interpretation and enrichment analysis.

> [!NOTE]
> The selected MCL granularity-4 solution balances strong cluster separation with high within-network density.

---

## 10. Network Analytics and Enrichment

Network-analysis outputs are located under:

```text
bone_EOA_analysis/network_reconstruction/
```

and include:

```text
network_analytics/
MCL_gran_4_cluster_enrichment/
```

Network analytics were performed using Cytoscape-exported network information together with Python/NetworkX-based analysis.

Calculated properties include:

- network density;
- node degree;
- clustering coefficient;
- betweenness centrality;
- closeness centrality;
- cluster density;
- cluster inertia; and
- silhouette-based cluster separation.

Node2Vec representations were additionally used for cluster-separation assessment based on the Cytoscape-derived cluster labels.

### Functional enrichment

Per-cluster enrichment analysis was performed on the selected MCL granularity-4 network.

The analyses include:

- Gene Ontology enrichment;
- Biological Process;
- Molecular Function;
- Cellular Component;
- protein-group enrichment; and
- CORUM-related enrichment analyses.

Where applicable, enrichment results were filtered using:

```text
Adjusted P-value < 0.05
```

The repository retains the corresponding enrichment tables, plots, and analysis code.

---

# Key Dataset Sizes

| Dataset | Rows | Role |
|---|---:|---|
| `bone_PPI_combs_raw.csv` | 877,501 | Raw OP PPI features |
| `bone_PPI_combs_processed.csv` | 877,501 | Processed OP PPI features |
| `bone_interactor_combs_raw.csv` | 494,515 | Raw interactor-combination features |
| `bone_interactor_combs_processed.csv` | 494,515 | Processed interactor-combination features |
| `bone_ppis_EOA_all_predictions_mapped.csv` | 877,501 | Complete mapped OP PPI predictions |
| `bone_interactor_combs_EOA_all_predictions_mapped.csv` | 494,515 | Complete mapped interactor predictions |
| `BOTH_FILTER_bone_ppis_EOA_positive_fully_mapped.csv` | 1,079 | Filtered OP PPIs |
| `BOTH_FILTER_bone_interactor_combs_EOA_positive_fully_mapped.csv` | 3,284 | Filtered OP-interactor PPIs |
| `final_bone_ppi_network_file.csv` | **4,363** | Final network |

---

# Dataset Columns

The final mapped interaction datasets contain the following principal fields:

| Column | Description |
|---|---|
| `uidA`, `uidB` | UniProt IDs of the interacting proteins |
| `Predicted Classes` | EOA predicted PPI class: `0` negative, `1` positive |
| `Probability Score` | Classifier prediction probability |
| `Regression Value` | EOA predicted interaction affinity |
| `mean_prob_aff` | Combined classifier-probability / predicted-affinity score |
| `stringdb_check` | Whether the pair occurs in the mapped STRING dataset |
| `score_norm` | Normalized STRING interaction score |
| `irefindex_check` | Whether the pair occurs in mapped iRefIndex |
| `uidA_irefindex` | Protein A UniProt ID as represented in iRefIndex |
| `uidB_irefindex` | Protein B UniProt ID as represented in iRefIndex |
| `method` | Experimental PPI detection method from iRefIndex |
| `Host_organism_taxid` | Host-organism information from iRefIndex |
| `numParticipants` | Number of proteins participating in the interaction |
| `GeneA`, `GeneB` | Gene names corresponding to proteins A and B |
| `OP_A`, `OP_B` | Whether protein A/B belongs to the OP list |
| `OP_check` | Distinguishes OP PPIs from expanded interactor PPIs |

---

# Repository Structure

```text
bone_interactome_analysis/
│
├── README.md
├── .gitignore
│
└── bone_EOA_analysis/
    │
    ├── Datasets/
    │   │
    │   ├── reference/
    │   │   ├── bone_proteins_04.xlsx
    │   │   ├── bone_comb.csv
    │   │   ├── uid_to_Gene_Names.csv
    │   │   ├── irefindex_v3.csv
    │   │   └── stringdb_ppis_curated.csv        [Google Drive]
    │   │
    │   ├── feature_calculated/                  [Google Drive]
    │   │   ├── bone_PPI_combs_raw.csv
    │   │   ├── bone_PPI_combs_processed.csv
    │   │   ├── bone_interactor_combs_raw.csv
    │   │   └── bone_interactor_combs_processed.csv
    │   │
    │   ├── intermediate/                        [Google Drive]
    │   │
    │   └── results/
    │       ├── bone_ppis_EOA_all_predictions_mapped.csv
    │       │                                           [Google Drive]
    │       ├── bone_interactor_combs_EOA_all_predictions_mapped.csv
    │       ├── BOTH_FILTER_bone_ppis_EOA_positive_fully_mapped.csv
    │       ├── BOTH_FILTER_bone_interactor_combs_EOA_positive_fully_mapped.csv
    │       └── final_bone_ppi_network_file.csv
    │
    ├── fc_EOA_prediction_codes/
    │   ├── prediction code
    │   └── models/
    │
    ├── mapping_codes/
    │   ├── bone_ppi_mappings.ipynb
    │   └── bone_interactor_combs_mappings.ipynb
    │
    └── network_reconstruction/
        │
        ├── bone_network.cys                       [Google Drive]
        ├── network figures
        │
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

> [!NOTE]
> Files marked `[Google Drive]` are intentionally excluded from Git because of their size. They should be restored to the indicated paths when reproducing the complete workflow.

---

# Large Files and Datasets

The complete large-file archive for this project is available here:

### [Download large datasets and Cytoscape project from Google Drive](PASTE_GOOGLE_DRIVE_FOLDER_LINK_HERE)

> [!WARNING]
> Large datasets are intentionally not tracked by Git.
>
> Do **not** commit these files directly to the repository. Download them from Google Drive and restore them to the paths shown below.

The external archive contains:

| Google Drive directory | Contents | Restore to |
|---|---|---|
| `feature_calculated/` | Raw and processed EOA feature matrices | `bone_EOA_analysis/Datasets/feature_calculated/` |
| `intermediate/` | Intermediate prediction and mapping datasets | `bone_EOA_analysis/Datasets/intermediate/` |
| `reference/` | Large external/reference datasets | `bone_EOA_analysis/Datasets/reference/` |
| `results/` | Large final EOA prediction datasets excluded from Git | `bone_EOA_analysis/Datasets/results/` |
| `cytoscape/` | Cytoscape project used for network reconstruction | `bone_EOA_analysis/network_reconstruction/` |

### Large files excluded from Git

Important excluded files include:

```text
Datasets/feature_calculated/
    bone_PPI_combs_raw.csv
    bone_PPI_combs_processed.csv
    bone_interactor_combs_raw.csv
    bone_interactor_combs_processed.csv

Datasets/reference/
    stringdb_ppis_curated.csv

Datasets/intermediate/
    [intermediate EOA prediction and mapping datasets]

Datasets/results/
    bone_ppis_EOA_all_predictions_mapped.csv

network_reconstruction/
    bone_network.cys
```

> [!TIP]
> After downloading the archive, preserve the directory structure shown above. The analysis code uses these directories to organize the different stages of the workflow.

---

# Reproducing the Analysis

A conceptual reproduction of Part 1 follows this order:

```text
44 osteoporosis-related proteins
              │
              ▼
20,421 reviewed human UniProt proteins
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
Probability > 0.517074
AND mean_prob_aff > 0.459257
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
Same dual filtering
              │
              ▼
3,284 filtered interactor PPIs
              │
              ▼
1,079 + 3,284
              │
              ▼
4,363 final PPIs
              │
              ▼
STRING / iRefIndex / Gene / OP mapping
              │
              ▼
Cytoscape network reconstruction
              │
              ▼
MCL / MCODE / Leiden / GLay comparison
              │
              ▼
MCL granularity 4
              │
              ▼
Network analytics + cluster enrichment
```

> [!IMPORTANT]
> Large datasets required for complete reproduction must first be downloaded from the [Google Drive archive](PASTE_GOOGLE_DRIVE_FOLDER_LINK_HERE) and restored to their corresponding directories.

---

# Acknowledgements

The present work has been developed as part of the **REGENERATION project**, funded by the European Union’s Horizon 2020 research and innovation program under the **Marie Sklodowska-Curie RISE (Grant Agreement No. 101131255)**.

This work was supported by the **Swiss State Secretariat for Education, Research and Innovation (SERI)** under contract number **23.0086**.
