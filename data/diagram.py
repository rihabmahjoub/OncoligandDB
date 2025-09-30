import pandas as pd
import matplotlib.pyplot as plt
from itertools import chain
from upsetplot import UpSet
import os

# --- Fichiers cancers solides (15) ---
solid_cancers = [
    ('breast.xlsx', 'Breast'),
    ('cervical.xlsx', 'Cervical'),
    ('bladder.xlsx', 'Bladder'),
    ('colorectal.xlsx', 'Colorectal'),
    ('esophageal.xlsx', 'Esophageal'),
    ('gastric.xlsx', 'Gastric'),
    ('head_and_neck.xlsx', 'Head and Neck'),
    ('hepatitis.xlsx', 'Hepatitis'),
    ('lung.xlsx', 'Lung'),
    ('melanoma.xlsx', 'Melanoma'),
    ('ovarian.xlsx', 'Ovarian'),
    ('pancreas.xlsx', 'Pancreas'),
    ('prostate.xlsx', 'Prostate'),
    ('renal.xlsx', 'Renal'),
    ('thyroid.xlsx', 'Thyroid')
]

# --- Fichiers cancers liquides (7) ---
liquid_cancers = [
    ('non_hodgkin_lymphoma.xlsx', 'Non-Hodgkin Lymphoma'),
    ('mds.xlsx', 'Myelodysplastic Syndromes (MDS)'),
    ('mm.xlsx', 'Multiple Myeloma (MM)'),
    ('cml.xlsx', 'Chronic Myeloid Leukemia (CML)'),
    ('cll.xlsx', 'Chronic Lymphocytic Leukemia (CLL)'),
    ('aml.xlsx', 'Acute Myeloid Leukemia (AML)'),
    ('all.xlsx', 'Acute Lymphoblastic Leukemia (ALL)')
]

# --- Fonction pour charger les ligands ---
def load_cancer_sets(cancer_files):
    cancer_sets = {}
    for file, cancer in cancer_files:
        if os.path.exists(file):
            df = pd.read_excel(file)
            # Normaliser les colonnes
            df.columns = df.columns.str.lower().str.replace(" ", "_")
            
            # Vérifier colonne du ligand
            if "ligand_name" in df.columns:
                ligands = set(df["ligand_name"].dropna().unique())
                cancer_sets[cancer] = ligands
            else:
                print(f"⚠️ Column 'Ligand Name' not found in {file}, skipping...")
        else:
            print(f"⚠️ File not found: {file}, skipping...")
    return cancer_sets

# --- Construire matrice binaire ---
def build_binary_df(cancer_sets):
    all_ligands = set(chain.from_iterable(cancer_sets.values()))
    records = []
    for ligand in all_ligands:
        row = {cancer: ligand in ligands for cancer, ligands in cancer_sets.items()}
        row["Ligand"] = ligand
        records.append(row)
    return pd.DataFrame(records).set_index(list(cancer_sets.keys()))

# --- UpSet plot ---
def plot_upset(cancer_sets, title):
    if not cancer_sets:
        print(f"❌ No data for {title}, skipping plot.")
        return
    df = build_binary_df(cancer_sets)
    upset = UpSet(df, subset_size='count', show_counts=True, sort_by="cardinality")
    upset.plot()
    plt.suptitle(title, fontsize=14, fontweight="bold")
    plt.show()

# --- Générer graphiques ---
solid_sets = load_cancer_sets(solid_cancers)
liquid_sets = load_cancer_sets(liquid_cancers)

if solid_sets:
    plot_upset(solid_sets, "Ligand Overlap Across Solid Cancers (15 types)")

if liquid_sets:
    plot_upset(liquid_sets, "Ligand Overlap Across Liquid Cancers (7 types)")

