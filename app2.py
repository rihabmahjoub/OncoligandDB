from flask import Flask, render_template, send_file, request, redirect, url_for 
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField
from wtforms.validators import DataRequired

import pandas as pd

import os

import openbabel

from rdkit import Chem

from rdkit.Chem import Draw



app = Flask(__name__)

app.config['SECRET_KEY'] = 'your_secret_key'



# Chemins des répertoires

data_dir = 'data/'

pdb_output_dir = "static/pdb"

image_output_dir = "static/images"



# Création des répertoires s'ils n'existent pas

os.makedirs(image_output_dir, exist_ok=True)

os.makedirs(pdb_output_dir, exist_ok=True)



# Fonction de lecture des données

def read_data(cancer_type):

    file_path = os.path.join(data_dir, f"{cancer_type}.xlsx")

    return pd.read_excel(file_path)



# Conversion SMILES vers PDB

def convert_smiles_to_pdb(smiles, output_file):

    obConversion = openbabel.OBConversion()

    obConversion.SetInAndOutFormats("smi", "pdb")

    mol = openbabel.OBMol()

    obConversion.ReadString(mol, smiles)

    obConversion.WriteFile(mol, output_file)



# Génération des images 2D

def generate_2d_images(df, output_dir=image_output_dir):

    for index, row in df.iterrows():

        smiles = row['SMILES']

        name = str(row['Product']).replace(" ", "_")

        mol = Chem.MolFromSmiles(smiles)

        if mol:

            path = os.path.join(output_dir, f"{name}.png")

            Draw.MolToFile(mol, path)

    return df



# Génération des liens de téléchargement

def generate_download_links(df, output_dir=pdb_output_dir):

    pdb_links = []

    for index, row in df.iterrows():

        smiles = row.get('SMILES')

        ligand_name = row.get('Product', f"ligand_{index}").replace(" ", "_")

        if pd.notna(smiles):

            pdb_filename = f"{ligand_name}.pdb"

            pdb_path = os.path.join(output_dir, pdb_filename)

            convert_smiles_to_pdb(smiles, pdb_path)

            pdb_links.append(f"/download_pdb?filename={pdb_filename}")

        else:

            pdb_links.append(None)

    df['pdb_link'] = pdb_links

    return df



# Formulaire SMILES

class SmilesForm(FlaskForm):

    smiles = StringField('SMILES', validators=[DataRequired()])

    ligand_name = StringField('Ligand Name', validators=[DataRequired()])

    submit = SubmitField('Convert to PDB')



# Routes

@app.route('/')

def index():

    return render_template('index.html')



@app.route('/solid_cancer')

def solid_cancer():

    return render_template('solid_cancer.html')



@app.route('/cancer_type/<cancer_type>')

def cancer_type(cancer_type):

    page = request.args.get('page', 1, type=int)

    per_page = 10



    df = read_data(cancer_type)

    df = generate_download_links(df)

    df = generate_2d_images(df)

    

    total = len(df)

    start = (page - 1) * per_page

    end = min(start + per_page, total)

    paginated_data = df.iloc[start:end].to_dict(orient='records')

    total_pages = (total + per_page - 1) // per_page



    form = SmilesForm()

    return render_template('cancer_type.html',

                           data=paginated_data,

                           cancer_type=cancer_type,

                           page=page,

                           total_pages=total_pages,

                           form=form)



@app.route('/convert_smiles_to_pdb', methods=['POST'])

def convert_smiles_to_pdb_route():

    form = SmilesForm()

    if form.validate_on_submit():

        smiles = form.smiles.data

        ligand_name = form.ligand_name.data.replace(" ", "_")

        output_file = f"{pdb_output_dir}/{ligand_name}.pdb"

        convert_smiles_to_pdb(smiles, output_file)

        return redirect(request.referrer)



@app.route('/download_pdb')

def download_pdb():

    filename = request.args.get('filename')

    return send_file(f"{pdb_output_dir}/{filename}", as_attachment=True)



if __name__ == '__main__':

    app.run(debug=True)

