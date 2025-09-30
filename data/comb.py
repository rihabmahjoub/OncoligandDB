import os
import pandas as pd

# Directory containing Excel files
data_dir = os.path.dirname(os.path.abspath(__file__))
print("Files in data_dir:", os.listdir(data_dir))  # debug check

# --- Solid cancers (15) ---
files = [
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
    ('thyroid.xlsx', 'Thyroid'),

    # --- Liquid cancers (7) ---
    ('non_hodgkin_lymphoma.xlsx', 'Non-Hodgkin Lymphoma'),
    ('mds.xlsx', 'Myelodysplastic Syndromes (MDS)'),
    ('mm.xlsx', 'Multiple Myeloma (MM)'),
    ('cml.xlsx', 'Chronic Myeloid Leukemia (CML)'),
    ('cll.xlsx', 'Chronic Lymphocytic Leukemia (CLL)'),
    ('aml.xlsx', 'Acute Myeloid Leukemia (AML)'),
    ('all.xlsx', 'Acute Lymphoblastic Leukemia (ALL)')
]

dfs = []

# Read each file and add Cancer_Type column
for file, cancer_type in files:
    file_path = os.path.join(data_dir, file)
    if not os.path.exists(file_path):
        print(f"⚠️ File not found: {file}, skipping...")
        continue
    
    df = pd.read_excel(file_path)
    df['Cancer_Type'] = cancer_type

    # Standardize column names
    df.columns = df.columns.str.lower().str.replace(' ', '_')

    # Handle product / ligand name
    if 'ligand_name' in df.columns:
        df.rename(columns={'ligand_name': 'product'}, inplace=True)
    elif 'name' in df.columns:
        df.rename(columns={'name': 'product'}, inplace=True)
    elif 'drug' in df.columns:
        df.rename(columns={'drug': 'product'}, inplace=True)
    elif 'product' not in df.columns:
        # fallback: create a placeholder product column
        df['product'] = [f"ligand_{i}" for i in range(len(df))]

    dfs.append(df)

# Combine DataFrames
if dfs:
    combined_df = pd.concat(dfs, ignore_index=True)

    # Standardize column order
    column_order = [
        'product', 'year', 'drugbank', 'chembl', 'indication', 'smiles', 'cancer_type'
    ]
    combined_df = combined_df.reindex(columns=column_order)

    out_path = os.path.join(data_dir, "combined_data.xlsx")
    combined_df.to_excel(out_path, index=False)

    # Summary by cancer type
    summary = combined_df.groupby("cancer_type").size()
    print("\n✅ Combined data saved to", out_path, f"({len(combined_df)} rows)")
    print("📊 Ligands per cancer type:\n", summary)
else:
    print("❌ No data combined (all files missing?)")

