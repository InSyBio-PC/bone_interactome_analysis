#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Protein UID Mapping Script

This script maps UniProt UIDs from an input dataset
(single column or PPI pairwise format) to cleaned UniProt/RefSeq/Gene info.
"""

import pandas as pd
import re
import os

# =========================== #
# USER CONFIGURATION
# =========================== #
single_uid_mode = 0  # <-- SET TO 1 if mapping single column of UIDs ('uid' column must exist),
                     #     0 if mapping PPI dataset ('uidA' & 'uidB' must exist)

input_file = (
    r"C:\Users\path\new_affinity_ds_2025_STRUCTURAL_PLUS_SEQUENTIAL_FEATURES.csv"
    if single_uid_mode == 0
    else r"C:\Users\harry\dataset.csv"
)

output_file = (
    r"C:\Users\path\dataset_mapped.csv"
    if single_uid_mode == 0
    else r"C:\Users\harry\dataset_mapped.csv"
)

mapping_file = r"C:\path\uniprotkb_HUMAN_REVIEWED_GENES_REFSEQ_SEQUENCES_2025_07_01.tsv"

# =========================== #
# LOAD AND CLEAN MAPPING DATA
# =========================== #

def pick_single_np(refseq_str):
    """Extract first NP_ RefSeq accession, prefer lowest -suffix number."""
    if pd.isna(refseq_str):
        return ''

    entries = [e.strip() for e in re.split(r'[;,]', str(refseq_str))]
    parsed = []
    for entry in entries:
        match = re.match(r'(NP_\d+\.\d+)(?:\s*\[(.*?)\])?', entry)
        if match:
            np_id = re.sub(r'\.\d+$', '', match.group(1))  # strip version
            annotation = match.group(2)
            parsed.append((np_id, annotation))

    if not parsed:
        return ''

    def get_suffix(annotation):
        if annotation and '-' in annotation:
            try:
                return int(annotation.split('-')[-1])
            except ValueError:
                return float('inf')
        return float('inf')  # no -x → least preferred

    parsed_sorted = sorted(parsed, key=lambda x: get_suffix(x[1]))
    return parsed_sorted[0][0] if parsed_sorted else ''


def pick_first_gene(gene_str):
    """Extract first gene name from UniProt 'Gene Names' field."""
    if pd.isna(gene_str):
        return ''
    return str(gene_str).strip().split()[0]


print("🔄 Loading UniProt mapping file...")
raw_df = pd.read_table(mapping_file, delimiter='\t')

# Clean RefSeq & Gene columns
raw_df['RefSeq_cleaned'] = raw_df['RefSeq'].apply(pick_single_np)
raw_df['Gene_cleaned'] = raw_df['Gene Names'].apply(pick_first_gene)

# Fallback: if gene missing, use UniProt Entry
raw_df['Gene_cleaned'] = raw_df.apply(
    lambda row: row['Entry'] if pd.isna(row['Gene_cleaned']) or row['Gene_cleaned'] == '' else row['Gene_cleaned'],
    axis=1
)

# Final mapping DataFrame
map_df = raw_df[['Entry', 'RefSeq_cleaned', 'Sequence', 'Gene_cleaned']].copy()
map_df.rename(
    columns={
        'Entry': 'uid',
        'RefSeq_cleaned': 'acession',
        'Sequence': 'seq',
        'Gene_cleaned': 'gene',
    },
    inplace=True,
)

print(f"✅ Mapping file loaded. {len(map_df)} entries available.")

# =========================== #
# APPLY MAPPING
# =========================== #

print("🔄 Loading user dataset...")
user_df = pd.read_csv(input_file)
print(f"✅ Input dataset loaded. Shape: {user_df.shape}")

# Strip whitespace from object columns
user_df = user_df.apply(lambda col: col.str.strip() if col.dtype == "object" else col)

if single_uid_mode == 0:
    # === Pairwise uidA / uidB mode ===
    if not {'uidA', 'uidB'}.issubset(user_df.columns):
        raise ValueError("In pairwise mode, input file must contain 'uidA' and 'uidB' columns.")

    mapped = user_df.merge(
        map_df.rename(columns={
            'uid': 'uidA',
            'acession': 'protein_accession_A',
            'seq': 'seq_A',
            'gene': 'GeneA',
        }),
        on='uidA',
        how='left'
    )

    mapped = mapped.merge(
        map_df.rename(columns={
            'uid': 'uidB',
            'acession': 'protein_accession_B',
            'seq': 'seq_B',
            'gene': 'GeneB',
        }),
        on='uidB',
        how='left'
    )

else:
    # === Single uid mode ===
    if 'uid' not in user_df.columns:
        raise ValueError("In single_uid_mode, input file must contain a 'uid' column.")

    mapped = user_df.merge(map_df, on='uid', how='left')

# Remove duplicates & suffix columns
mapped = mapped.drop_duplicates()
mapped = mapped.loc[:, ~mapped.columns.str.endswith(('_x', '_y'))]

# Save result
os.makedirs(os.path.dirname(output_file), exist_ok=True)
mapped.to_csv(output_file, index=False)

print(f"✅ Mapping complete. Output saved to:\n{output_file}")
print(f"📊 Final dataset shape: {mapped.shape}")
