from flask import Flask, render_template, send_file, request, redirect, url_for
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField
from wtforms.validators import DataRequired
import pandas as pd
import os
import openbabel

app = Flask(__name__)
app.config['SECRET_KEY'] = 'your_secret_key'

# Path to local database (Excel files)
data_dir = 'data/'

# Function to read Excel data
def read_data(cancer_type):
    file_path = os.path.join(data_dir, f"{cancer_type}.xlsx")
    df = pd.read_excel(file_path)
    return df

# Function to convert SMILES to PDB
def convert_smiles_to_pdb(smiles, output_file):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "pdb")

    mol = openbabel.OBMol()
    obConversion.ReadString(mol, smiles)
    obConversion.WriteFile(mol, output_file)
    return output_file

# Generate PDB download links
def generate_download_links(df, output_dir="static/pdb"):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    pdb_links = []
    for index, row in df.iterrows():
        smiles = row.get('SMILES')
        ligand_name = row.get('Product', f"ligand_{index}")
        if pd.notna(smiles):
            pdb_filename = f"{output_dir}/{ligand_name}.pdb"
            convert_smiles_to_pdb(smiles, pdb_filename)
            pdb_links.append(f"/download_pdb?filename={ligand_name}.pdb")
        else:
            pdb_links.append(None)

    df['pdb_link'] = pdb_links
    return df

# Flask-WTF form
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
    df = read_data(cancer_type)
    df = generate_download_links(df)
    data = df.to_dict(orient='records')
    form = SmilesForm()
    return render_template('cancer_type.html', data=data, cancer_type=cancer_type, form=form)

@app.route('/convert_smiles_to_pdb', methods=['POST'])
def convert_smiles_to_pdb_route():
    form = SmilesForm()
    if form.validate_on_submit():
        smiles = form.smiles.data
        ligand_name = form.ligand_name.data
        output_file = f"static/pdb/{ligand_name}.pdb"
        convert_smiles_to_pdb(smiles, output_file)
    return redirect(request.referrer)

@app.route('/download_pdb')
def download_pdb():
    filename = request.args.get('filename')
    return send_file(f"static/pdb/{filename}", as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)

