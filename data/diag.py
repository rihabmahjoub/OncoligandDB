import pandas as pd
import matplotlib.pyplot as plt
from itertools import chain
from upsetplot import UpSet
import seaborn as sns

# --- Apply seaborn style ---
sns.set(style="whitegrid")

# --- Load combined dataset ---
df_all = pd.read_excel("combined_data.xlsx")
df_all = df_all.dropna(subset=["product"])  # drop missing ligands

# --- Define solid and liquid cancers ---
solid_cancers = [
    "Breast", "Cervical", "Bladder", "Colorec", "Esophag",
    "Gastric", "H&N", "Hepatitis", "Lung", "Melan",
    "Ovarian", "Pancreas", "Prostate", "Renal", "Thyroid"
]

liquid_cancers = [
    "Non-Hodgkin Lymphoma", "Myelodysplastic Syndromes (MDS)",
    "Multiple Myeloma (MM)", "Chronic Myeloid Leukemia (CML)",
    "Chronic Lymphocytic Leukemia (CLL)",
    "Acute Myeloid Leukemia (AML)", "Acute Lymphoblastic Leukemia (ALL)"
]

def build_upset(df, cancer_list, title, figsize=(14,7)):
    """Build and plot an UpSet diagram for selected cancers"""
    # --- Build ligand sets ---
    cancer_sets = {
        cancer: set(df.loc[df['cancer_type'] == cancer, 'product'].unique())
        for cancer in cancer_list
    }

    # --- Create binary matrix ---
    all_ligands = set(chain.from_iterable(cancer_sets.values()))
    records = []
    for ligand in all_ligands:
        row = {cancer: ligand in ligands for cancer, ligands in cancer_sets.items()}
        row["Ligand"] = ligand
        records.append(row)

    df_matrix = pd.DataFrame(records).set_index(cancer_list)

    # --- Plot ---
    plt.figure(figsize=figsize)
    upset = UpSet(
        df_matrix,
        subset_size='count',
        show_counts=True,
        orientation='horizontal',
        sort_by='cardinality'
    )
    upset.plot()

    # Title
    plt.suptitle(title, fontsize=16, fontweight='bold', color='#333333')
    plt.gcf().patch.set_facecolor('#f5f5f5')
    plt.tight_layout()
    plt.show()


# --- Plot solid cancers ---
build_upset(df_all, solid_cancers, "Ligand Overlap Across 15 Solid Cancer Types")

# --- Plot liquid cancers ---
build_upset(df_all, liquid_cancers, "Ligand Overlap Across 7 Liquid Cancer Types")

