import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Charger le fichier combiné
df = pd.read_excel("combined_data.xlsx")

# Créer la figure et les axes avec plus de largeur pour les noms longs
fig, ax = plt.subplots(figsize=(12, 10))

# --- Plot: Nombre de ligands par type de cancer ---
if "cancer_type" in df.columns:
    sns.countplot(
        y="cancer_type",
        data=df,
        order=df['cancer_type'].value_counts().index,
        palette="viridis",
        ax=ax
    )
    ax.set_title("Number of Ligands per Cancer Type", fontsize=14, fontweight='bold')
    ax.set_xlabel("Number of Ligands")
    ax.set_ylabel("Cancer Type")
    ax.grid(axis='x', linestyle='--', alpha=0.7)

    # Ajuster l'affichage des noms longs
    plt.tight_layout()  # évite que les labels soient coupés
else:
    ax.text(0.5, 0.5, "⚠️ 'cancer_type' column missing", ha="center", va="center")

# ✅ Sauvegarde en PNG
plt.savefig("ligands_per_cancer_type_barplot.png", dpi=300, bbox_inches="tight")
plt.show()

