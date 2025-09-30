import pandas as pd
from shiny import App, render, ui
import plotly.express as px

# --- Load dataset ---
df = pd.read_excel("data/combined_data.xlsx")

# --- Clean column names ---
df.columns = df.columns.str.strip().str.lower().str.replace(" ", "_")

# --- Add PubChem Image column ---
def smiles_to_pubchem_img(smiles):
    if pd.isna(smiles) or smiles == "not found":
        return ""
    return f'<img src="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/PNG" width="120">'

df["2D_image"] = df["smiles"].apply(smiles_to_pubchem_img)

# --- Shiny UI ---
app_ui = ui.page_fluid(
    ui.h2("OncoligandDB: Dashboard of Cancer Ligands"),
    ui.navset_tab(
        ui.nav_panel("Ligand Table",
            ui.output_ui("ligand_table")
        ),
        ui.nav_panel("Plots",
            ui.layout_columns(
                ui.card(ui.output_plot("year_plot")),
                ui.card(ui.output_plot("smiles_plot")),
                ui.card(ui.output_plot("cancer_plot"))
            )
        )
    )
)

# --- Shiny Server ---
def server(input, output, session):

    @output
    @render.ui
    def ligand_table():
        # Render table with images
        return ui.HTML(df.to_html(escape=False, index=False))

    @output
    @render.plot
    def year_plot():
        if "year" not in df.columns:
            return
        fig = px.histogram(df.dropna(subset=["year"]), x="year", nbins=30,
                           title="Distribution of Drug Approval Years")
        return fig

    @output
    @render.plot
    def smiles_plot():
        if "smiles" not in df.columns:
            return
        smiles_data = df[df["smiles"].notna() & (df["smiles"] != "not found")].copy()
        smiles_data["length"] = smiles_data["smiles"].apply(len)
        fig = px.histogram(smiles_data, x="length", nbins=30,
                           title="Distribution of SMILES Lengths")
        return fig

    @output
    @render.plot
    def cancer_plot():
        if "cancer_type" not in df.columns:
            return
        fig = px.histogram(df, y="cancer_type", title="Number of Ligands per Cancer Type")
        return fig

# --- Run App ---
app = App(app_ui, server)

