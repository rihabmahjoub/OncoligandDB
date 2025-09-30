import pandas as pd
import glob
import os

# 📂 Folder containing your 22 Excel files
data_dir = "/home/rihab/Téléchargements/oncoligandDB/data"
files = glob.glob(os.path.join(data_dir, "*.xlsx"))

# List to store all data
all_data = []

for file in files:
    try:
        df = pd.read_excel(file)

        # Normalize columns: lowercase + replace spaces with underscores
        df.columns = df.columns.str.lower().str.strip().str.replace(" ", "_")

        # Check if 'ligand_name' exists
        if "ligand_name" in df.columns:
            # Select relevant columns if they exist
            relevant_cols = ["ligand_name", "year", "drugbank", "chembl", "indication", "smiles", "2d_image", "pdb_link"]
            existing_cols = [col for col in relevant_cols if col in df.columns]
            df_subset = df[existing_cols].dropna(subset=["ligand_name"])
            all_data.append(df_subset)
        else:
            print(f"⚠️ Column 'ligand_name' not found in {os.path.basename(file)}, skipped...")

    except Exception as e:
        print(f"❌ Error with {os.path.basename(file)}: {e}")

# Concatenate all data
if all_data:
    full_df = pd.concat(all_data, ignore_index=True)

    # Remove duplicate ligands (keep first occurrence)
    full_df_unique = full_df.drop_duplicates(subset=["ligand_name"])

    # Save to CSV
    output_csv = os.path.join(data_dir, "unique_ligands_full.csv")
    full_df_unique.to_csv(output_csv, index=False)
    print(f"✅ Total unique ligands: {len(full_df_unique)}")
    print(f"💾 Saved full table to: {output_csv}")
else:
    print("❌ No valid data found in any Excel files.")

