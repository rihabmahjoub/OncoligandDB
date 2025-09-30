import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Charger le fichier combiné
df = pd.read_excel("combined_data.xlsx")

# Nettoyer les colonnes (tout en minuscule + enlever espaces)
df.columns = df.columns.str.lower().str.replace(" ", "_")

# 🔹 Nettoyage de la colonne 'year'
if "year" in df.columns:
    df['year'] = pd.to_numeric(df['year'], errors='coerce')  # convertit en numérique, invalide -> NaN
    df = df.dropna(subset=['year'])  # supprimer les NaN dans year
    df['year'] = df['year'].astype(int)  # convertir en entier

# --- Préparer la figure globale ---
fig, axes = plt.subplots(2, 2, figsize=(16, 12))
plt.subplots_adjust(hspace=0.4, wspace=0.3)  # espace entre subplots

# --- Plot 1: Distribution des années d’approbation ---
if "year" in df.columns:
    sns.histplot(data=df, x='year', bins=30, kde=True, color='blue', ax=axes[0, 0])
    axes[0, 0].set_title('Distribution of Drug Approval Years', fontsize=14, fontweight='bold')
    axes[0, 0].set_xlabel('Year')
    axes[0, 0].set_ylabel('Number of Drugs')
    axes[0, 0].grid(axis='y', linestyle='--', alpha=0.7)
else:
    axes[0, 0].text(0.5, 0.5, "⚠️ 'year' column missing", ha="center", va="center")

# --- Plot 2: Distribution des longueurs SMILES ---
if "smiles" in df.columns:
    smiles_data = df[df['smiles'].notna() & (df['smiles'] != 'not found')].copy()
    smiles_data['smiles_length'] = smiles_data['smiles'].apply(len)
    sns.histplot(data=smiles_data, x='smiles_length', bins=30, kde=True, color='salmon', ax=axes[0, 1])
    axes[0, 1].set_title('Distribution of SMILES String Lengths', fontsize=14, fontweight='bold')
    axes[0, 1].set_xlabel('SMILES Length (Characters)')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].grid(axis='y', linestyle='--', alpha=0.7)
else:
    axes[0, 1].text(0.5, 0.5, "⚠️ 'smiles' column missing", ha="center", va="center")

# --- Plot 3: Nombre de ligands par type de cancer ---
if "cancer_type" in df.columns:
    sns.countplot(y="cancer_type", data=df, order=df['cancer_type'].value_counts().index, palette="viridis", ax=axes[1, 0])
    axes[1, 0].set_title("Number of Ligands per Cancer Type", fontsize=14, fontweight='bold')
    axes[1, 0].set_xlabel("Number of Ligands")
    axes[1, 0].set_ylabel("Cancer Type")
    axes[1, 0].grid(axis='x', linestyle='--', alpha=0.7)
else:
    axes[1, 0].text(0.5, 0.5, "⚠️ 'cancer_type' column missing", ha="center", va="center")

# --- Plot 4: Top 15 DrugBank IDs les plus fréquents ---
if "drugbank" in df.columns:
    top_drugs = df['drugbank'].value_counts().head(15)
    sns.barplot(x=top_drugs.values, y=top_drugs.index, palette="mako", ax=axes[1, 1])
    axes[1, 1].set_title("Top 15 Most Frequent DrugBank IDs", fontsize=14, fontweight='bold')
    axes[1, 1].set_xlabel("Count")
    axes[1, 1].set_ylabel("DrugBank ID")
    axes[1, 1].grid(axis='x', linestyle='--', alpha=0.7)
else:
    axes[1, 1].text(0.5, 0.5, "⚠️ 'drugbank' column missing", ha="center", va="center")

# ✅ Sauvegarde en PNG
plt.savefig("oncoliganddb_summary.png", dpi=300, bbox_inches="tight")
plt.show()

